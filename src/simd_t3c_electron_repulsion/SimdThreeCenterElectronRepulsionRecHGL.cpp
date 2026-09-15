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


#include "SimdThreeCenterElectronRepulsionRecHGL.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSL.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
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
#include "SimdTransferHG.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hgl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hgl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 218851, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1683 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 218851, 175183, 10535, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 17,
                                                             ncols, fj, 6, fq);

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

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 77, 0, 3, 8, 9,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 83, 0, 3, 9, 10,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 89, 0, 3, 10, 11,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 95, 0, 3, 11, 12,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 101, 0, 3, 12, 13,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 107, 0, 3, 13, 14,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 14, 15,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 119, 0, 3, 15, 16,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 125, 0, 3, 16, 17,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 131, 0, 3, 17, 18,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 137, 0, 3, 18, 19,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 19, 20,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 149, 0, 3, 20, 21,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 155, 0, 3, 21, 22,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 161, 0, 3, 22, 23,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 167, 0, 3, 23, 24,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 26, 29,
                                                                       77, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 29, 32,
                                                                       83, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 32, 35,
                                                                       89, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 35, 38,
                                                                       95, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 38, 41,
                                                                       101, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 223, 0, 3, 41, 44,
                                                                       107, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 233, 0, 3, 44, 47,
                                                                       113, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 243, 0, 3, 47, 50,
                                                                       119, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 253, 0, 3, 50, 53,
                                                                       125, 131, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 263, 0, 3, 53, 56,
                                                                       131, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 56, 59,
                                                                       137, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 283, 0, 3, 59, 62,
                                                                       143, 149, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 293, 0, 3, 62, 65,
                                                                       149, 155, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 303, 0, 3, 65, 68,
                                                                       155, 161, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 313, 0, 3, 68, 71,
                                                                       161, 167, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 323, 0, 3, 77, 83,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 83, 89,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 353, 0, 3, 89, 95,
                                                                       193, 203, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 95,
                                                                       101, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 383, 0, 3, 101,
                                                                       107, 213, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 107,
                                                                       113, 223, 233, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 413, 0, 3, 113,
                                                                       119, 233, 243, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 119,
                                                                       125, 243, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 443, 0, 3, 125,
                                                                       131, 253, 263, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 131,
                                                                       137, 263, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 473, 0, 3, 137,
                                                                       143, 273, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 488, 0, 3, 143,
                                                                       149, 283, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 503, 0, 3, 149,
                                                                       155, 293, 303, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 518, 0, 3, 155,
                                                                       161, 303, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 533, 0, 3, 173,
                                                                       183, 323, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 554, 0, 3, 183,
                                                                       193, 338, 353, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 575, 0, 3, 193,
                                                                       203, 353, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 596, 0, 3, 203,
                                                                       213, 368, 383, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 617, 0, 3, 213,
                                                                       223, 383, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 638, 0, 3, 223,
                                                                       233, 398, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 659, 0, 3, 233,
                                                                       243, 413, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 680, 0, 3, 243,
                                                                       253, 428, 443, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 701, 0, 3, 253,
                                                                       263, 443, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 722, 0, 3, 263,
                                                                       273, 458, 473, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 743, 0, 3, 273,
                                                                       283, 473, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 764, 0, 3, 283,
                                                                       293, 488, 503, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 785, 0, 3, 293,
                                                                       303, 503, 518, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 806, 0, 3, 323,
                                                                       338, 533, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 834, 0, 3, 338,
                                                                       353, 554, 575, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 862, 0, 3, 353,
                                                                       368, 575, 596, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 890, 0, 3, 368,
                                                                       383, 596, 617, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 918, 0, 3, 383,
                                                                       398, 617, 638, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 946, 0, 3, 398,
                                                                       413, 638, 659, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 974, 0, 3, 413,
                                                                       428, 659, 680, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1002, 0, 3, 428,
                                                                       443, 680, 701, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 443,
                                                                       458, 701, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1058, 0, 3, 458,
                                                                       473, 722, 743, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1086, 0, 3, 473,
                                                                       488, 743, 764, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1114, 0, 3, 488,
                                                                       503, 764, 785, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1142, 0, 3, 533,
                                                                       554, 806, 834, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1178, 0, 3, 554,
                                                                       575, 834, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1214, 0, 3, 575,
                                                                       596, 862, 890, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1250, 0, 3, 596,
                                                                       617, 890, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1286, 0, 3, 617,
                                                                       638, 918, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1322, 0, 3, 638,
                                                                       659, 946, 974, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1358, 0, 3, 659,
                                                                       680, 974, 1002, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1394, 0, 3, 680,
                                                                       701, 1002, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1430, 0, 3, 701,
                                                                       722, 1030, 1058, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1466, 0, 3, 722,
                                                                       743, 1058, 1086, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1502, 0, 3, 743,
                                                                       764, 1086, 1114, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1538, 0, 3, 806,
                                                                       834, 1142, 1178, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1583, 0, 3, 834,
                                                                       862, 1178, 1214, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1628, 0, 3, 862,
                                                                       890, 1214, 1250, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1673, 0, 3, 890,
                                                                       918, 1250, 1286, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1718, 0, 3, 918,
                                                                       946, 1286, 1322, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1763, 0, 3, 946,
                                                                       974, 1322, 1358, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1808, 0, 3, 974,
                                                                       1002, 1358, 1394, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1853, 0, 3, 1002,
                                                                       1030, 1394, 1430, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1898, 0, 3, 1030,
                                                                       1058, 1430, 1466, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1943, 0, 3, 1058,
                                                                       1086, 1466, 1502, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1988, 0, 3, 1142,
                                                                       1178, 1538, 1583, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2043, 0, 3, 1178,
                                                                       1214, 1583, 1628, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2098, 0, 3, 1214,
                                                                       1250, 1628, 1673, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2153, 0, 3, 1250,
                                                                       1286, 1673, 1718, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2208, 0, 3, 1286,
                                                                       1322, 1718, 1763, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2263, 0, 3, 1322,
                                                                       1358, 1763, 1808, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2318, 0, 3, 1358,
                                                                       1394, 1808, 1853, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2373, 0, 3, 1394,
                                                                       1430, 1853, 1898, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 2428, 0, 3, 1430,
                                                                       1466, 1898, 1943, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2483, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2486, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2489, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2492, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2495, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2498, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2501, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2504, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2507, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2510, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2513, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2516, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2519, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2522, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2525, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2528, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2531, 3, 10, 32,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2540, 3, 11, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2549, 3, 12, 38,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2558, 3, 13, 41,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2567, 3, 14, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2576, 3, 15, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2585, 3, 16, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2594, 3, 17, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2603, 3, 18, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2612, 3, 19, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2621, 3, 20, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2630, 3, 21, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2639, 3, 22, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2648, 3, 23, 71,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2657, 3, 24, 74,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2666, 3, 32, 89,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2684, 3, 35, 95,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2702, 3, 38, 101,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2720, 3, 41, 107,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2738, 3, 44, 113,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2756, 3, 47, 119,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2774, 3, 50, 125,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2792, 3, 53, 131,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2810, 3, 56, 137,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2828, 3, 59, 143,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2846, 3, 62, 149,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2864, 3, 65, 155,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2882, 3, 68, 161,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2900, 3, 71, 167,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2918, 3, 89, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2948, 3, 95, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2978, 3, 101, 213,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3008, 3, 107, 223,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3038, 3, 113, 233,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3068, 3, 119, 243,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3098, 3, 125, 253,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3128, 3, 131, 263,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3158, 3, 137, 273,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3188, 3, 143, 283,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3218, 3, 149, 293,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3248, 3, 155, 303,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3278, 3, 161, 313,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3308, 3, 193, 353,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3353, 3, 203, 368,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3398, 3, 213, 383,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3443, 3, 223, 398,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3488, 3, 233, 413,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3533, 3, 243, 428,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3578, 3, 253, 443,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3623, 3, 263, 458,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3668, 3, 273, 473,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3713, 3, 283, 488,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3758, 3, 293, 503,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3803, 3, 303, 518,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3848, 3, 353, 575,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3911, 3, 368, 596,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3974, 3, 383, 617,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4037, 3, 398, 638,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4100, 3, 413, 659,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4163, 3, 428, 680,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4226, 3, 443, 701,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4289, 3, 458, 722,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4352, 3, 473, 743,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4415, 3, 488, 764,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4478, 3, 503, 785,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4541, 3, 575, 862,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4625, 3, 596, 890,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4709, 3, 617, 918,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4793, 3, 638, 946,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4877, 3, 659, 974,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4961, 3, 680,
                                                                       1002, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5045, 3, 701,
                                                                       1030, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5129, 3, 722,
                                                                       1058, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5213, 3, 743,
                                                                       1086, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5297, 3, 764,
                                                                       1114, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5381, 3, 862,
                                                                       1214, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5489, 3, 890,
                                                                       1250, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5597, 3, 918,
                                                                       1286, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5705, 3, 946,
                                                                       1322, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5813, 3, 974,
                                                                       1358, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5921, 3, 1002,
                                                                       1394, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6029, 3, 1030,
                                                                       1430, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6137, 3, 1058,
                                                                       1466, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6245, 3, 1086,
                                                                       1502, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6353, 3, 1214,
                                                                       1628, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6488, 3, 1250,
                                                                       1673, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6623, 3, 1286,
                                                                       1718, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6758, 3, 1322,
                                                                       1763, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6893, 3, 1358,
                                                                       1808, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7028, 3, 1394,
                                                                       1853, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7163, 3, 1430,
                                                                       1898, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7298, 3, 1466,
                                                                       1943, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7433, 3, 1628,
                                                                       2098, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7598, 3, 1673,
                                                                       2153, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7763, 3, 1718,
                                                                       2208, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7928, 3, 1763,
                                                                       2263, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8093, 3, 1808,
                                                                       2318, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8258, 3, 1853,
                                                                       2373, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8423, 3, 1898,
                                                                       2428, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8588, 3, 8, 9,
                                                                       2483, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8594, 3, 9, 10,
                                                                       2486, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8600, 3, 10, 11,
                                                                       2489, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8606, 3, 11, 12,
                                                                       2492, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8612, 3, 12, 13,
                                                                       2495, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8618, 3, 13, 14,
                                                                       2498, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8624, 3, 14, 15,
                                                                       2501, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8630, 3, 15, 16,
                                                                       2504, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8636, 3, 16, 17,
                                                                       2507, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8642, 3, 17, 18,
                                                                       2510, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8648, 3, 18, 19,
                                                                       2513, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8654, 3, 19, 20,
                                                                       2516, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8660, 3, 20, 21,
                                                                       2519, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8666, 3, 21, 22,
                                                                       2522, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8672, 3, 22, 23,
                                                                       2525, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8678, 3, 23, 24,
                                                                       2528, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8684, 0, 3, 8588,
                                                                       2483, 8594, 2531, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8702, 0, 3, 8594,
                                                                       2486, 8600, 2540, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8720, 0, 3, 8600,
                                                                       2489, 8606, 2549, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8738, 0, 3, 8606,
                                                                       2492, 8612, 2558, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8756, 0, 3, 8612,
                                                                       2495, 8618, 2567, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8774, 0, 3, 8618,
                                                                       2498, 8624, 2576, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8792, 0, 3, 8624,
                                                                       2501, 8630, 2585, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8810, 0, 3, 8630,
                                                                       2504, 8636, 2594, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8828, 0, 3, 8636,
                                                                       2507, 8642, 2603, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8846, 0, 3, 8642,
                                                                       2510, 8648, 2612, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8864, 0, 3, 8648,
                                                                       2513, 8654, 2621, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8882, 0, 3, 8654,
                                                                       2516, 8660, 2630, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8900, 0, 3, 8660,
                                                                       2519, 8666, 2639, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8918, 0, 3, 8666,
                                                                       2522, 8672, 2648, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8936, 0, 3, 8672,
                                                                       2525, 8678, 2657, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8954, 0, 3, 8684,
                                                                       2531, 8702, 77, 83, 2666,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8990, 0, 3, 8702,
                                                                       2540, 8720, 83, 89, 2684,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9026, 0, 3, 8720,
                                                                       2549, 8738, 89, 95, 2702,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9062, 0, 3, 8738,
                                                                       2558, 8756, 95, 101, 2720,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9098, 0, 3, 8756,
                                                                       2567, 8774, 101, 107,
                                                                       2738, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9134, 0, 3, 8774,
                                                                       2576, 8792, 107, 113,
                                                                       2756, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9170, 0, 3, 8792,
                                                                       2585, 8810, 113, 119,
                                                                       2774, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9206, 0, 3, 8810,
                                                                       2594, 8828, 119, 125,
                                                                       2792, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9242, 0, 3, 8828,
                                                                       2603, 8846, 125, 131,
                                                                       2810, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9278, 0, 3, 8846,
                                                                       2612, 8864, 131, 137,
                                                                       2828, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9314, 0, 3, 8864,
                                                                       2621, 8882, 137, 143,
                                                                       2846, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9350, 0, 3, 8882,
                                                                       2630, 8900, 143, 149,
                                                                       2864, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9386, 0, 3, 8900,
                                                                       2639, 8918, 149, 155,
                                                                       2882, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9422, 0, 3, 8918,
                                                                       2648, 8936, 155, 161,
                                                                       2900, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9458, 0, 3, 8954,
                                                                       2666, 8990, 173, 183,
                                                                       2918, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9518, 0, 3, 8990,
                                                                       2684, 9026, 183, 193,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9578, 0, 3, 9026,
                                                                       2702, 9062, 193, 203,
                                                                       2978, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9638, 0, 3, 9062,
                                                                       2720, 9098, 203, 213,
                                                                       3008, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9698, 0, 3, 9098,
                                                                       2738, 9134, 213, 223,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9758, 0, 3, 9134,
                                                                       2756, 9170, 223, 233,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9818, 0, 3, 9170,
                                                                       2774, 9206, 233, 243,
                                                                       3098, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9878, 0, 3, 9206,
                                                                       2792, 9242, 243, 253,
                                                                       3128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9938, 0, 3, 9242,
                                                                       2810, 9278, 253, 263,
                                                                       3158, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9998, 0, 3, 9278,
                                                                       2828, 9314, 263, 273,
                                                                       3188, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10058, 0, 3, 9314,
                                                                       2846, 9350, 273, 283,
                                                                       3218, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10118, 0, 3, 9350,
                                                                       2864, 9386, 283, 293,
                                                                       3248, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10178, 0, 3, 9386,
                                                                       2882, 9422, 293, 303,
                                                                       3278, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10238, 0, 3, 9458,
                                                                       2918, 9518, 323, 338,
                                                                       3308, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10328, 0, 3, 9518,
                                                                       2948, 9578, 338, 353,
                                                                       3353, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10418, 0, 3, 9578,
                                                                       2978, 9638, 353, 368,
                                                                       3398, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10508, 0, 3, 9638,
                                                                       3008, 9698, 368, 383,
                                                                       3443, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10598, 0, 3, 9698,
                                                                       3038, 9758, 383, 398,
                                                                       3488, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10688, 0, 3, 9758,
                                                                       3068, 9818, 398, 413,
                                                                       3533, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10778, 0, 3, 9818,
                                                                       3098, 9878, 413, 428,
                                                                       3578, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10868, 0, 3, 9878,
                                                                       3128, 9938, 428, 443,
                                                                       3623, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10958, 0, 3, 9938,
                                                                       3158, 9998, 443, 458,
                                                                       3668, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11048, 0, 3, 9998,
                                                                       3188, 10058, 458, 473,
                                                                       3713, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11138, 0, 3,
                                                                       10058, 3218, 10118, 473,
                                                                       488, 3758, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11228, 0, 3,
                                                                       10118, 3248, 10178, 488,
                                                                       503, 3803, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11318, 0, 3,
                                                                       10238, 3308, 10328, 533,
                                                                       554, 3848, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11444, 0, 3,
                                                                       10328, 3353, 10418, 554,
                                                                       575, 3911, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11570, 0, 3,
                                                                       10418, 3398, 10508, 575,
                                                                       596, 3974, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11696, 0, 3,
                                                                       10508, 3443, 10598, 596,
                                                                       617, 4037, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11822, 0, 3,
                                                                       10598, 3488, 10688, 617,
                                                                       638, 4100, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11948, 0, 3,
                                                                       10688, 3533, 10778, 638,
                                                                       659, 4163, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12074, 0, 3,
                                                                       10778, 3578, 10868, 659,
                                                                       680, 4226, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12200, 0, 3,
                                                                       10868, 3623, 10958, 680,
                                                                       701, 4289, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12326, 0, 3,
                                                                       10958, 3668, 11048, 701,
                                                                       722, 4352, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12452, 0, 3,
                                                                       11048, 3713, 11138, 722,
                                                                       743, 4415, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12578, 0, 3,
                                                                       11138, 3758, 11228, 743,
                                                                       764, 4478, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12704, 0, 3,
                                                                       11318, 3848, 11444, 806,
                                                                       834, 4541, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12872, 0, 3,
                                                                       11444, 3911, 11570, 834,
                                                                       862, 4625, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13040, 0, 3,
                                                                       11570, 3974, 11696, 862,
                                                                       890, 4709, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13208, 0, 3,
                                                                       11696, 4037, 11822, 890,
                                                                       918, 4793, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13376, 0, 3,
                                                                       11822, 4100, 11948, 918,
                                                                       946, 4877, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13544, 0, 3,
                                                                       11948, 4163, 12074, 946,
                                                                       974, 4961, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13712, 0, 3,
                                                                       12074, 4226, 12200, 974,
                                                                       1002, 5045, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13880, 0, 3,
                                                                       12200, 4289, 12326, 1002,
                                                                       1030, 5129, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14048, 0, 3,
                                                                       12326, 4352, 12452, 1030,
                                                                       1058, 5213, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14216, 0, 3,
                                                                       12452, 4415, 12578, 1058,
                                                                       1086, 5297, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14384, 0, 3,
                                                                       12704, 4541, 12872, 1142,
                                                                       1178, 5381, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14600, 0, 3,
                                                                       12872, 4625, 13040, 1178,
                                                                       1214, 5489, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14816, 0, 3,
                                                                       13040, 4709, 13208, 1214,
                                                                       1250, 5597, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15032, 0, 3,
                                                                       13208, 4793, 13376, 1250,
                                                                       1286, 5705, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15248, 0, 3,
                                                                       13376, 4877, 13544, 1286,
                                                                       1322, 5813, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15464, 0, 3,
                                                                       13544, 4961, 13712, 1322,
                                                                       1358, 5921, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15680, 0, 3,
                                                                       13712, 5045, 13880, 1358,
                                                                       1394, 6029, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15896, 0, 3,
                                                                       13880, 5129, 14048, 1394,
                                                                       1430, 6137, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16112, 0, 3,
                                                                       14048, 5213, 14216, 1430,
                                                                       1466, 6245, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16328, 0, 3,
                                                                       14384, 5381, 14600, 1538,
                                                                       1583, 6353, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16598, 0, 3,
                                                                       14600, 5489, 14816, 1583,
                                                                       1628, 6488, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16868, 0, 3,
                                                                       14816, 5597, 15032, 1628,
                                                                       1673, 6623, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17138, 0, 3,
                                                                       15032, 5705, 15248, 1673,
                                                                       1718, 6758, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17408, 0, 3,
                                                                       15248, 5813, 15464, 1718,
                                                                       1763, 6893, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17678, 0, 3,
                                                                       15464, 5921, 15680, 1763,
                                                                       1808, 7028, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17948, 0, 3,
                                                                       15680, 6029, 15896, 1808,
                                                                       1853, 7163, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 18218, 0, 3,
                                                                       15896, 6137, 16112, 1853,
                                                                       1898, 7298, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 18488, 0, 3,
                                                                       16328, 6353, 16598, 1988,
                                                                       2043, 7433, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 18818, 0, 3,
                                                                       16598, 6488, 16868, 2043,
                                                                       2098, 7598, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19148, 0, 3,
                                                                       16868, 6623, 17138, 2098,
                                                                       2153, 7763, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19478, 0, 3,
                                                                       17138, 6758, 17408, 2153,
                                                                       2208, 7928, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 19808, 0, 3,
                                                                       17408, 6893, 17678, 2208,
                                                                       2263, 8093, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20138, 0, 3,
                                                                       17678, 7028, 17948, 2263,
                                                                       2318, 8258, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 20468, 0, 3,
                                                                       17948, 7163, 18218, 2318,
                                                                       2373, 8423, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20798, 3, 2483,
                                                                       2486, 8600, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20808, 3, 2486,
                                                                       2489, 8606, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20818, 3, 2489,
                                                                       2492, 8612, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20828, 3, 2492,
                                                                       2495, 8618, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20838, 3, 2495,
                                                                       2498, 8624, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20848, 3, 2498,
                                                                       2501, 8630, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20858, 3, 2501,
                                                                       2504, 8636, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20868, 3, 2504,
                                                                       2507, 8642, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20878, 3, 2507,
                                                                       2510, 8648, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20888, 3, 2510,
                                                                       2513, 8654, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20898, 3, 2513,
                                                                       2516, 8660, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20908, 3, 2516,
                                                                       2519, 8666, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20918, 3, 2519,
                                                                       2522, 8672, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 20928, 3, 2522,
                                                                       2525, 8678, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 20938, 0, 3,
                                                                       20798, 8600, 20808, 2531,
                                                                       2540, 8720, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 20968, 0, 3,
                                                                       20808, 8606, 20818, 2540,
                                                                       2549, 8738, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 20998, 0, 3,
                                                                       20818, 8612, 20828, 2549,
                                                                       2558, 8756, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21028, 0, 3,
                                                                       20828, 8618, 20838, 2558,
                                                                       2567, 8774, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21058, 0, 3,
                                                                       20838, 8624, 20848, 2567,
                                                                       2576, 8792, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21088, 0, 3,
                                                                       20848, 8630, 20858, 2576,
                                                                       2585, 8810, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21118, 0, 3,
                                                                       20858, 8636, 20868, 2585,
                                                                       2594, 8828, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21148, 0, 3,
                                                                       20868, 8642, 20878, 2594,
                                                                       2603, 8846, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21178, 0, 3,
                                                                       20878, 8648, 20888, 2603,
                                                                       2612, 8864, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21208, 0, 3,
                                                                       20888, 8654, 20898, 2612,
                                                                       2621, 8882, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21238, 0, 3,
                                                                       20898, 8660, 20908, 2621,
                                                                       2630, 8900, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21268, 0, 3,
                                                                       20908, 8666, 20918, 2630,
                                                                       2639, 8918, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21298, 0, 3,
                                                                       20918, 8672, 20928, 2639,
                                                                       2648, 8936, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 21328, 0, 3,
                                                                       20938, 8720, 20968, 2666,
                                                                       2684, 9026, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 21388, 0, 3,
                                                                       20968, 8738, 20998, 2684,
                                                                       2702, 9062, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 21448, 0, 3,
                                                                       20998, 8756, 21028, 2702,
                                                                       2720, 9098, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 21508, 0, 3,
                                                                       21028, 8774, 21058, 2720,
                                                                       2738, 9134, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 21568, 0, 3,
                                                                       21058, 8792, 21088, 2738,
                                                                       2756, 9170, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 21628, 0, 3,
                                                                       21088, 8810, 21118, 2756,
                                                                       2774, 9206, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 21688, 0, 3,
                                                                       21118, 8828, 21148, 2774,
                                                                       2792, 9242, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 21748, 0, 3,
                                                                       21148, 8846, 21178, 2792,
                                                                       2810, 9278, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 21808, 0, 3,
                                                                       21178, 8864, 21208, 2810,
                                                                       2828, 9314, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 21868, 0, 3,
                                                                       21208, 8882, 21238, 2828,
                                                                       2846, 9350, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 21928, 0, 3,
                                                                       21238, 8900, 21268, 2846,
                                                                       2864, 9386, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 21988, 0, 3,
                                                                       21268, 8918, 21298, 2864,
                                                                       2882, 9422, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22048, 0, 3,
                                                                       21328, 9026, 21388, 2918,
                                                                       2948, 9578, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22148, 0, 3,
                                                                       21388, 9062, 21448, 2948,
                                                                       2978, 9638, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22248, 0, 3,
                                                                       21448, 9098, 21508, 2978,
                                                                       3008, 9698, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22348, 0, 3,
                                                                       21508, 9134, 21568, 3008,
                                                                       3038, 9758, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22448, 0, 3,
                                                                       21568, 9170, 21628, 3038,
                                                                       3068, 9818, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22548, 0, 3,
                                                                       21628, 9206, 21688, 3068,
                                                                       3098, 9878, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22648, 0, 3,
                                                                       21688, 9242, 21748, 3098,
                                                                       3128, 9938, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22748, 0, 3,
                                                                       21748, 9278, 21808, 3128,
                                                                       3158, 9998, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22848, 0, 3,
                                                                       21808, 9314, 21868, 3158,
                                                                       3188, 10058, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22948, 0, 3,
                                                                       21868, 9350, 21928, 3188,
                                                                       3218, 10118, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23048, 0, 3,
                                                                       21928, 9386, 21988, 3218,
                                                                       3248, 10178, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23148, 0, 3,
                                                                       22048, 9578, 22148, 3308,
                                                                       3353, 10418, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23298, 0, 3,
                                                                       22148, 9638, 22248, 3353,
                                                                       3398, 10508, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23448, 0, 3,
                                                                       22248, 9698, 22348, 3398,
                                                                       3443, 10598, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23598, 0, 3,
                                                                       22348, 9758, 22448, 3443,
                                                                       3488, 10688, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23748, 0, 3,
                                                                       22448, 9818, 22548, 3488,
                                                                       3533, 10778, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23898, 0, 3,
                                                                       22548, 9878, 22648, 3533,
                                                                       3578, 10868, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24048, 0, 3,
                                                                       22648, 9938, 22748, 3578,
                                                                       3623, 10958, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24198, 0, 3,
                                                                       22748, 9998, 22848, 3623,
                                                                       3668, 11048, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24348, 0, 3,
                                                                       22848, 10058, 22948, 3668,
                                                                       3713, 11138, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24498, 0, 3,
                                                                       22948, 10118, 23048, 3713,
                                                                       3758, 11228, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 24648, 0, 3,
                                                                       23148, 10418, 23298, 3848,
                                                                       3911, 11570, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 24858, 0, 3,
                                                                       23298, 10508, 23448, 3911,
                                                                       3974, 11696, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25068, 0, 3,
                                                                       23448, 10598, 23598, 3974,
                                                                       4037, 11822, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25278, 0, 3,
                                                                       23598, 10688, 23748, 4037,
                                                                       4100, 11948, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25488, 0, 3,
                                                                       23748, 10778, 23898, 4100,
                                                                       4163, 12074, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25698, 0, 3,
                                                                       23898, 10868, 24048, 4163,
                                                                       4226, 12200, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25908, 0, 3,
                                                                       24048, 10958, 24198, 4226,
                                                                       4289, 12326, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26118, 0, 3,
                                                                       24198, 11048, 24348, 4289,
                                                                       4352, 12452, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26328, 0, 3,
                                                                       24348, 11138, 24498, 4352,
                                                                       4415, 12578, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 26538, 0, 3,
                                                                       24648, 11570, 24858, 4541,
                                                                       4625, 13040, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 26818, 0, 3,
                                                                       24858, 11696, 25068, 4625,
                                                                       4709, 13208, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27098, 0, 3,
                                                                       25068, 11822, 25278, 4709,
                                                                       4793, 13376, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27378, 0, 3,
                                                                       25278, 11948, 25488, 4793,
                                                                       4877, 13544, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27658, 0, 3,
                                                                       25488, 12074, 25698, 4877,
                                                                       4961, 13712, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27938, 0, 3,
                                                                       25698, 12200, 25908, 4961,
                                                                       5045, 13880, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28218, 0, 3,
                                                                       25908, 12326, 26118, 5045,
                                                                       5129, 14048, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28498, 0, 3,
                                                                       26118, 12452, 26328, 5129,
                                                                       5213, 14216, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 28778, 0, 3,
                                                                       26538, 13040, 26818, 5381,
                                                                       5489, 14816, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 29138, 0, 3,
                                                                       26818, 13208, 27098, 5489,
                                                                       5597, 15032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 29498, 0, 3,
                                                                       27098, 13376, 27378, 5597,
                                                                       5705, 15248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 29858, 0, 3,
                                                                       27378, 13544, 27658, 5705,
                                                                       5813, 15464, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 30218, 0, 3,
                                                                       27658, 13712, 27938, 5813,
                                                                       5921, 15680, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 30578, 0, 3,
                                                                       27938, 13880, 28218, 5921,
                                                                       6029, 15896, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 30938, 0, 3,
                                                                       28218, 14048, 28498, 6029,
                                                                       6137, 16112, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 31298, 0, 3,
                                                                       28778, 14816, 29138, 6353,
                                                                       6488, 16868, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 31748, 0, 3,
                                                                       29138, 15032, 29498, 6488,
                                                                       6623, 17138, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 32198, 0, 3,
                                                                       29498, 15248, 29858, 6623,
                                                                       6758, 17408, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 32648, 0, 3,
                                                                       29858, 15464, 30218, 6758,
                                                                       6893, 17678, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 33098, 0, 3,
                                                                       30218, 15680, 30578, 6893,
                                                                       7028, 17948, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 33548, 0, 3,
                                                                       30578, 15896, 30938, 7028,
                                                                       7163, 18218, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 33998, 0, 3,
                                                                       31298, 16868, 31748, 7433,
                                                                       7598, 19148, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 34548, 0, 3,
                                                                       31748, 17138, 32198, 7598,
                                                                       7763, 19478, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 35098, 0, 3,
                                                                       32198, 17408, 32648, 7763,
                                                                       7928, 19808, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 35648, 0, 3,
                                                                       32648, 17678, 33098, 7928,
                                                                       8093, 20138, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 36198, 0, 3,
                                                                       33098, 17948, 33548, 8093,
                                                                       8258, 20468, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36748, 3, 8588,
                                                                       8594, 20798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36763, 3, 8594,
                                                                       8600, 20808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36778, 3, 8600,
                                                                       8606, 20818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36793, 3, 8606,
                                                                       8612, 20828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36808, 3, 8612,
                                                                       8618, 20838, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36823, 3, 8618,
                                                                       8624, 20848, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36838, 3, 8624,
                                                                       8630, 20858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36853, 3, 8630,
                                                                       8636, 20868, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36868, 3, 8636,
                                                                       8642, 20878, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36883, 3, 8642,
                                                                       8648, 20888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36898, 3, 8648,
                                                                       8654, 20898, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36913, 3, 8654,
                                                                       8660, 20908, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36928, 3, 8660,
                                                                       8666, 20918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 36943, 3, 8666,
                                                                       8672, 20928, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 36958, 0, 3,
                                                                       36748, 20798, 36763, 8684,
                                                                       8702, 20938, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 37003, 0, 3,
                                                                       36763, 20808, 36778, 8702,
                                                                       8720, 20968, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 37048, 0, 3,
                                                                       36778, 20818, 36793, 8720,
                                                                       8738, 20998, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 37093, 0, 3,
                                                                       36793, 20828, 36808, 8738,
                                                                       8756, 21028, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 37138, 0, 3,
                                                                       36808, 20838, 36823, 8756,
                                                                       8774, 21058, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 37183, 0, 3,
                                                                       36823, 20848, 36838, 8774,
                                                                       8792, 21088, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 37228, 0, 3,
                                                                       36838, 20858, 36853, 8792,
                                                                       8810, 21118, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 37273, 0, 3,
                                                                       36853, 20868, 36868, 8810,
                                                                       8828, 21148, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 37318, 0, 3,
                                                                       36868, 20878, 36883, 8828,
                                                                       8846, 21178, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 37363, 0, 3,
                                                                       36883, 20888, 36898, 8846,
                                                                       8864, 21208, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 37408, 0, 3,
                                                                       36898, 20898, 36913, 8864,
                                                                       8882, 21238, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 37453, 0, 3,
                                                                       36913, 20908, 36928, 8882,
                                                                       8900, 21268, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 37498, 0, 3,
                                                                       36928, 20918, 36943, 8900,
                                                                       8918, 21298, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37543, 0, 3,
                                                                       36958, 20938, 37003, 8954,
                                                                       8990, 21328, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37633, 0, 3,
                                                                       37003, 20968, 37048, 8990,
                                                                       9026, 21388, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37723, 0, 3,
                                                                       37048, 20998, 37093, 9026,
                                                                       9062, 21448, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37813, 0, 3,
                                                                       37093, 21028, 37138, 9062,
                                                                       9098, 21508, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37903, 0, 3,
                                                                       37138, 21058, 37183, 9098,
                                                                       9134, 21568, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 37993, 0, 3,
                                                                       37183, 21088, 37228, 9134,
                                                                       9170, 21628, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 38083, 0, 3,
                                                                       37228, 21118, 37273, 9170,
                                                                       9206, 21688, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 38173, 0, 3,
                                                                       37273, 21148, 37318, 9206,
                                                                       9242, 21748, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 38263, 0, 3,
                                                                       37318, 21178, 37363, 9242,
                                                                       9278, 21808, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 38353, 0, 3,
                                                                       37363, 21208, 37408, 9278,
                                                                       9314, 21868, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 38443, 0, 3,
                                                                       37408, 21238, 37453, 9314,
                                                                       9350, 21928, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 38533, 0, 3,
                                                                       37453, 21268, 37498, 9350,
                                                                       9386, 21988, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38623, 0, 3,
                                                                       37543, 21328, 37633, 9458,
                                                                       9518, 22048, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38773, 0, 3,
                                                                       37633, 21388, 37723, 9518,
                                                                       9578, 22148, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38923, 0, 3,
                                                                       37723, 21448, 37813, 9578,
                                                                       9638, 22248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 39073, 0, 3,
                                                                       37813, 21508, 37903, 9638,
                                                                       9698, 22348, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 39223, 0, 3,
                                                                       37903, 21568, 37993, 9698,
                                                                       9758, 22448, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 39373, 0, 3,
                                                                       37993, 21628, 38083, 9758,
                                                                       9818, 22548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 39523, 0, 3,
                                                                       38083, 21688, 38173, 9818,
                                                                       9878, 22648, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 39673, 0, 3,
                                                                       38173, 21748, 38263, 9878,
                                                                       9938, 22748, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 39823, 0, 3,
                                                                       38263, 21808, 38353, 9938,
                                                                       9998, 22848, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 39973, 0, 3,
                                                                       38353, 21868, 38443, 9998,
                                                                       10058, 22948, ncols,
                                                                       gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 40123, 0, 3,
                                                                       38443, 21928, 38533,
                                                                       10058, 10118, 23048,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 40273, 0, 3,
                                                                       38623, 22048, 38773,
                                                                       10238, 10328, 23148,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 40498, 0, 3,
                                                                       38773, 22148, 38923,
                                                                       10328, 10418, 23298,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 40723, 0, 3,
                                                                       38923, 22248, 39073,
                                                                       10418, 10508, 23448,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 40948, 0, 3,
                                                                       39073, 22348, 39223,
                                                                       10508, 10598, 23598,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 41173, 0, 3,
                                                                       39223, 22448, 39373,
                                                                       10598, 10688, 23748,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 41398, 0, 3,
                                                                       39373, 22548, 39523,
                                                                       10688, 10778, 23898,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 41623, 0, 3,
                                                                       39523, 22648, 39673,
                                                                       10778, 10868, 24048,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 41848, 0, 3,
                                                                       39673, 22748, 39823,
                                                                       10868, 10958, 24198,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 42073, 0, 3,
                                                                       39823, 22848, 39973,
                                                                       10958, 11048, 24348,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 42298, 0, 3,
                                                                       39973, 22948, 40123,
                                                                       11048, 11138, 24498,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 42523, 0, 3,
                                                                       40273, 23148, 40498,
                                                                       11318, 11444, 24648,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 42838, 0, 3,
                                                                       40498, 23298, 40723,
                                                                       11444, 11570, 24858,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 43153, 0, 3,
                                                                       40723, 23448, 40948,
                                                                       11570, 11696, 25068,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 43468, 0, 3,
                                                                       40948, 23598, 41173,
                                                                       11696, 11822, 25278,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 43783, 0, 3,
                                                                       41173, 23748, 41398,
                                                                       11822, 11948, 25488,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 44098, 0, 3,
                                                                       41398, 23898, 41623,
                                                                       11948, 12074, 25698,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 44413, 0, 3,
                                                                       41623, 24048, 41848,
                                                                       12074, 12200, 25908,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 44728, 0, 3,
                                                                       41848, 24198, 42073,
                                                                       12200, 12326, 26118,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 45043, 0, 3,
                                                                       42073, 24348, 42298,
                                                                       12326, 12452, 26328,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 45358, 0, 3,
                                                                       42523, 24648, 42838,
                                                                       12704, 12872, 26538,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 45778, 0, 3,
                                                                       42838, 24858, 43153,
                                                                       12872, 13040, 26818,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 46198, 0, 3,
                                                                       43153, 25068, 43468,
                                                                       13040, 13208, 27098,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 46618, 0, 3,
                                                                       43468, 25278, 43783,
                                                                       13208, 13376, 27378,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 47038, 0, 3,
                                                                       43783, 25488, 44098,
                                                                       13376, 13544, 27658,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 47458, 0, 3,
                                                                       44098, 25698, 44413,
                                                                       13544, 13712, 27938,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 47878, 0, 3,
                                                                       44413, 25908, 44728,
                                                                       13712, 13880, 28218,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 48298, 0, 3,
                                                                       44728, 26118, 45043,
                                                                       13880, 14048, 28498,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 48718, 0, 3,
                                                                       45358, 26538, 45778,
                                                                       14384, 14600, 28778,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 49258, 0, 3,
                                                                       45778, 26818, 46198,
                                                                       14600, 14816, 29138,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 49798, 0, 3,
                                                                       46198, 27098, 46618,
                                                                       14816, 15032, 29498,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 50338, 0, 3,
                                                                       46618, 27378, 47038,
                                                                       15032, 15248, 29858,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 50878, 0, 3,
                                                                       47038, 27658, 47458,
                                                                       15248, 15464, 30218,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 51418, 0, 3,
                                                                       47458, 27938, 47878,
                                                                       15464, 15680, 30578,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 51958, 0, 3,
                                                                       47878, 28218, 48298,
                                                                       15680, 15896, 30938,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 52498, 0, 3,
                                                                       48718, 28778, 49258,
                                                                       16328, 16598, 31298,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 53173, 0, 3,
                                                                       49258, 29138, 49798,
                                                                       16598, 16868, 31748,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 53848, 0, 3,
                                                                       49798, 29498, 50338,
                                                                       16868, 17138, 32198,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 54523, 0, 3,
                                                                       50338, 29858, 50878,
                                                                       17138, 17408, 32648,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 55198, 0, 3,
                                                                       50878, 30218, 51418,
                                                                       17408, 17678, 33098,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 55873, 0, 3,
                                                                       51418, 30578, 51958,
                                                                       17678, 17948, 33548,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 56548, 0, 3,
                                                                       52498, 31298, 53173,
                                                                       18488, 18818, 33998,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 57373, 0, 3,
                                                                       53173, 31748, 53848,
                                                                       18818, 19148, 34548,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 58198, 0, 3,
                                                                       53848, 32198, 54523,
                                                                       19148, 19478, 35098,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 59023, 0, 3,
                                                                       54523, 32648, 55198,
                                                                       19478, 19808, 35648,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 59848, 0, 3,
                                                                       55198, 33098, 55873,
                                                                       19808, 20138, 36198,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60673, 3, 20798,
                                                                       20808, 36778, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60694, 3, 20808,
                                                                       20818, 36793, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60715, 3, 20818,
                                                                       20828, 36808, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60736, 3, 20828,
                                                                       20838, 36823, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60757, 3, 20838,
                                                                       20848, 36838, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60778, 3, 20848,
                                                                       20858, 36853, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60799, 3, 20858,
                                                                       20868, 36868, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60820, 3, 20868,
                                                                       20878, 36883, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60841, 3, 20878,
                                                                       20888, 36898, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60862, 3, 20888,
                                                                       20898, 36913, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60883, 3, 20898,
                                                                       20908, 36928, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 60904, 3, 20908,
                                                                       20918, 36943, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 60925, 0, 3,
                                                                       60673, 36778, 60694,
                                                                       20938, 20968, 37048,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 60988, 0, 3,
                                                                       60694, 36793, 60715,
                                                                       20968, 20998, 37093,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 61051, 0, 3,
                                                                       60715, 36808, 60736,
                                                                       20998, 21028, 37138,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 61114, 0, 3,
                                                                       60736, 36823, 60757,
                                                                       21028, 21058, 37183,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 61177, 0, 3,
                                                                       60757, 36838, 60778,
                                                                       21058, 21088, 37228,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 61240, 0, 3,
                                                                       60778, 36853, 60799,
                                                                       21088, 21118, 37273,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 61303, 0, 3,
                                                                       60799, 36868, 60820,
                                                                       21118, 21148, 37318,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 61366, 0, 3,
                                                                       60820, 36883, 60841,
                                                                       21148, 21178, 37363,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 61429, 0, 3,
                                                                       60841, 36898, 60862,
                                                                       21178, 21208, 37408,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 61492, 0, 3,
                                                                       60862, 36913, 60883,
                                                                       21208, 21238, 37453,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 61555, 0, 3,
                                                                       60883, 36928, 60904,
                                                                       21238, 21268, 37498,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 61618, 0, 3,
                                                                       60925, 37048, 60988,
                                                                       21328, 21388, 37723,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 61744, 0, 3,
                                                                       60988, 37093, 61051,
                                                                       21388, 21448, 37813,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 61870, 0, 3,
                                                                       61051, 37138, 61114,
                                                                       21448, 21508, 37903,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 61996, 0, 3,
                                                                       61114, 37183, 61177,
                                                                       21508, 21568, 37993,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 62122, 0, 3,
                                                                       61177, 37228, 61240,
                                                                       21568, 21628, 38083,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 62248, 0, 3,
                                                                       61240, 37273, 61303,
                                                                       21628, 21688, 38173,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 62374, 0, 3,
                                                                       61303, 37318, 61366,
                                                                       21688, 21748, 38263,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 62500, 0, 3,
                                                                       61366, 37363, 61429,
                                                                       21748, 21808, 38353,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 62626, 0, 3,
                                                                       61429, 37408, 61492,
                                                                       21808, 21868, 38443,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 62752, 0, 3,
                                                                       61492, 37453, 61555,
                                                                       21868, 21928, 38533,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 62878, 0, 3,
                                                                       61618, 37723, 61744,
                                                                       22048, 22148, 38923,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 63088, 0, 3,
                                                                       61744, 37813, 61870,
                                                                       22148, 22248, 39073,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 63298, 0, 3,
                                                                       61870, 37903, 61996,
                                                                       22248, 22348, 39223,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 63508, 0, 3,
                                                                       61996, 37993, 62122,
                                                                       22348, 22448, 39373,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 63718, 0, 3,
                                                                       62122, 38083, 62248,
                                                                       22448, 22548, 39523,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 63928, 0, 3,
                                                                       62248, 38173, 62374,
                                                                       22548, 22648, 39673,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 64138, 0, 3,
                                                                       62374, 38263, 62500,
                                                                       22648, 22748, 39823,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 64348, 0, 3,
                                                                       62500, 38353, 62626,
                                                                       22748, 22848, 39973,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 64558, 0, 3,
                                                                       62626, 38443, 62752,
                                                                       22848, 22948, 40123,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 64768, 0, 3,
                                                                       62878, 38923, 63088,
                                                                       23148, 23298, 40723,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 65083, 0, 3,
                                                                       63088, 39073, 63298,
                                                                       23298, 23448, 40948,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 65398, 0, 3,
                                                                       63298, 39223, 63508,
                                                                       23448, 23598, 41173,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 65713, 0, 3,
                                                                       63508, 39373, 63718,
                                                                       23598, 23748, 41398,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 66028, 0, 3,
                                                                       63718, 39523, 63928,
                                                                       23748, 23898, 41623,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 66343, 0, 3,
                                                                       63928, 39673, 64138,
                                                                       23898, 24048, 41848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 66658, 0, 3,
                                                                       64138, 39823, 64348,
                                                                       24048, 24198, 42073,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 66973, 0, 3,
                                                                       64348, 39973, 64558,
                                                                       24198, 24348, 42298,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 67288, 0, 3,
                                                                       64768, 40723, 65083,
                                                                       24648, 24858, 43153,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 67729, 0, 3,
                                                                       65083, 40948, 65398,
                                                                       24858, 25068, 43468,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 68170, 0, 3,
                                                                       65398, 41173, 65713,
                                                                       25068, 25278, 43783,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 68611, 0, 3,
                                                                       65713, 41398, 66028,
                                                                       25278, 25488, 44098,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 69052, 0, 3,
                                                                       66028, 41623, 66343,
                                                                       25488, 25698, 44413,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 69493, 0, 3,
                                                                       66343, 41848, 66658,
                                                                       25698, 25908, 44728,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 69934, 0, 3,
                                                                       66658, 42073, 66973,
                                                                       25908, 26118, 45043,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 70375, 0, 3,
                                                                       67288, 43153, 67729,
                                                                       26538, 26818, 46198,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 70963, 0, 3,
                                                                       67729, 43468, 68170,
                                                                       26818, 27098, 46618,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 71551, 0, 3,
                                                                       68170, 43783, 68611,
                                                                       27098, 27378, 47038,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 72139, 0, 3,
                                                                       68611, 44098, 69052,
                                                                       27378, 27658, 47458,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 72727, 0, 3,
                                                                       69052, 44413, 69493,
                                                                       27658, 27938, 47878,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 73315, 0, 3,
                                                                       69493, 44728, 69934,
                                                                       27938, 28218, 48298,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 73903, 0, 3,
                                                                       70375, 46198, 70963,
                                                                       28778, 29138, 49798,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 74659, 0, 3,
                                                                       70963, 46618, 71551,
                                                                       29138, 29498, 50338,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 75415, 0, 3,
                                                                       71551, 47038, 72139,
                                                                       29498, 29858, 50878,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 76171, 0, 3,
                                                                       72139, 47458, 72727,
                                                                       29858, 30218, 51418,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 76927, 0, 3,
                                                                       72727, 47878, 73315,
                                                                       30218, 30578, 51958,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 77683, 0, 3,
                                                                       73903, 49798, 74659,
                                                                       31298, 31748, 53848,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 78628, 0, 3,
                                                                       74659, 50338, 75415,
                                                                       31748, 32198, 54523,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 79573, 0, 3,
                                                                       75415, 50878, 76171,
                                                                       32198, 32648, 55198,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 80518, 0, 3,
                                                                       76171, 51418, 76927,
                                                                       32648, 33098, 55873,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 81463, 0, 3,
                                                                       77683, 53848, 78628,
                                                                       33998, 34548, 58198,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 82618, 0, 3,
                                                                       78628, 54523, 79573,
                                                                       34548, 35098, 59023,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 83773, 0, 3,
                                                                       79573, 55198, 80518,
                                                                       35098, 35648, 59848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 84928, 3, 36748,
                                                                       36763, 60673, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 84956, 3, 36763,
                                                                       36778, 60694, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 84984, 3, 36778,
                                                                       36793, 60715, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85012, 3, 36793,
                                                                       36808, 60736, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85040, 3, 36808,
                                                                       36823, 60757, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85068, 3, 36823,
                                                                       36838, 60778, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85096, 3, 36838,
                                                                       36853, 60799, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85124, 3, 36853,
                                                                       36868, 60820, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85152, 3, 36868,
                                                                       36883, 60841, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85180, 3, 36883,
                                                                       36898, 60862, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85208, 3, 36898,
                                                                       36913, 60883, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 85236, 3, 36913,
                                                                       36928, 60904, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 85264, 0, 3,
                                                                       84928, 60673, 84956,
                                                                       36958, 37003, 60925,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 85348, 0, 3,
                                                                       84956, 60694, 84984,
                                                                       37003, 37048, 60988,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 85432, 0, 3,
                                                                       84984, 60715, 85012,
                                                                       37048, 37093, 61051,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 85516, 0, 3,
                                                                       85012, 60736, 85040,
                                                                       37093, 37138, 61114,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 85600, 0, 3,
                                                                       85040, 60757, 85068,
                                                                       37138, 37183, 61177,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 85684, 0, 3,
                                                                       85068, 60778, 85096,
                                                                       37183, 37228, 61240,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 85768, 0, 3,
                                                                       85096, 60799, 85124,
                                                                       37228, 37273, 61303,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 85852, 0, 3,
                                                                       85124, 60820, 85152,
                                                                       37273, 37318, 61366,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 85936, 0, 3,
                                                                       85152, 60841, 85180,
                                                                       37318, 37363, 61429,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 86020, 0, 3,
                                                                       85180, 60862, 85208,
                                                                       37363, 37408, 61492,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 86104, 0, 3,
                                                                       85208, 60883, 85236,
                                                                       37408, 37453, 61555,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 86188, 0, 3,
                                                                       85264, 60925, 85348,
                                                                       37543, 37633, 61618,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 86356, 0, 3,
                                                                       85348, 60988, 85432,
                                                                       37633, 37723, 61744,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 86524, 0, 3,
                                                                       85432, 61051, 85516,
                                                                       37723, 37813, 61870,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 86692, 0, 3,
                                                                       85516, 61114, 85600,
                                                                       37813, 37903, 61996,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 86860, 0, 3,
                                                                       85600, 61177, 85684,
                                                                       37903, 37993, 62122,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 87028, 0, 3,
                                                                       85684, 61240, 85768,
                                                                       37993, 38083, 62248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 87196, 0, 3,
                                                                       85768, 61303, 85852,
                                                                       38083, 38173, 62374,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 87364, 0, 3,
                                                                       85852, 61366, 85936,
                                                                       38173, 38263, 62500,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 87532, 0, 3,
                                                                       85936, 61429, 86020,
                                                                       38263, 38353, 62626,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 87700, 0, 3,
                                                                       86020, 61492, 86104,
                                                                       38353, 38443, 62752,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 87868, 0, 3,
                                                                       86188, 61618, 86356,
                                                                       38623, 38773, 62878,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 88148, 0, 3,
                                                                       86356, 61744, 86524,
                                                                       38773, 38923, 63088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 88428, 0, 3,
                                                                       86524, 61870, 86692,
                                                                       38923, 39073, 63298,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 88708, 0, 3,
                                                                       86692, 61996, 86860,
                                                                       39073, 39223, 63508,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 88988, 0, 3,
                                                                       86860, 62122, 87028,
                                                                       39223, 39373, 63718,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 89268, 0, 3,
                                                                       87028, 62248, 87196,
                                                                       39373, 39523, 63928,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 89548, 0, 3,
                                                                       87196, 62374, 87364,
                                                                       39523, 39673, 64138,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 89828, 0, 3,
                                                                       87364, 62500, 87532,
                                                                       39673, 39823, 64348,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 90108, 0, 3,
                                                                       87532, 62626, 87700,
                                                                       39823, 39973, 64558,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 90388, 0, 3,
                                                                       87868, 62878, 88148,
                                                                       40273, 40498, 64768,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 90808, 0, 3,
                                                                       88148, 63088, 88428,
                                                                       40498, 40723, 65083,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 91228, 0, 3,
                                                                       88428, 63298, 88708,
                                                                       40723, 40948, 65398,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 91648, 0, 3,
                                                                       88708, 63508, 88988,
                                                                       40948, 41173, 65713,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 92068, 0, 3,
                                                                       88988, 63718, 89268,
                                                                       41173, 41398, 66028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 92488, 0, 3,
                                                                       89268, 63928, 89548,
                                                                       41398, 41623, 66343,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 92908, 0, 3,
                                                                       89548, 64138, 89828,
                                                                       41623, 41848, 66658,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 93328, 0, 3,
                                                                       89828, 64348, 90108,
                                                                       41848, 42073, 66973,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 93748, 0, 3,
                                                                       90388, 64768, 90808,
                                                                       42523, 42838, 67288,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 94336, 0, 3,
                                                                       90808, 65083, 91228,
                                                                       42838, 43153, 67729,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 94924, 0, 3,
                                                                       91228, 65398, 91648,
                                                                       43153, 43468, 68170,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 95512, 0, 3,
                                                                       91648, 65713, 92068,
                                                                       43468, 43783, 68611,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 96100, 0, 3,
                                                                       92068, 66028, 92488,
                                                                       43783, 44098, 69052,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 96688, 0, 3,
                                                                       92488, 66343, 92908,
                                                                       44098, 44413, 69493,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 97276, 0, 3,
                                                                       92908, 66658, 93328,
                                                                       44413, 44728, 69934,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 97864, 0, 3,
                                                                       93748, 67288, 94336,
                                                                       45358, 45778, 70375,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 98648, 0, 3,
                                                                       94336, 67729, 94924,
                                                                       45778, 46198, 70963,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 99432, 0, 3,
                                                                       94924, 68170, 95512,
                                                                       46198, 46618, 71551,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 100216, 0, 3,
                                                                       95512, 68611, 96100,
                                                                       46618, 47038, 72139,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 101000, 0, 3,
                                                                       96100, 69052, 96688,
                                                                       47038, 47458, 72727,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 101784, 0, 3,
                                                                       96688, 69493, 97276,
                                                                       47458, 47878, 73315,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 102568, 0, 3,
                                                                       97864, 70375, 98648,
                                                                       48718, 49258, 73903,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 103576, 0, 3,
                                                                       98648, 70963, 99432,
                                                                       49258, 49798, 74659,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 104584, 0, 3,
                                                                       99432, 71551, 100216,
                                                                       49798, 50338, 75415,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 105592, 0, 3,
                                                                       100216, 72139, 101000,
                                                                       50338, 50878, 76171,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 106600, 0, 3,
                                                                       101000, 72727, 101784,
                                                                       50878, 51418, 76927,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 107608, 0, 3,
                                                                       102568, 73903, 103576,
                                                                       52498, 53173, 77683,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 108868, 0, 3,
                                                                       103576, 74659, 104584,
                                                                       53173, 53848, 78628,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 110128, 0, 3,
                                                                       104584, 75415, 105592,
                                                                       53848, 54523, 79573,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 111388, 0, 3,
                                                                       105592, 76171, 106600,
                                                                       54523, 55198, 80518,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 112648, 0, 3,
                                                                       107608, 77683, 108868,
                                                                       56548, 57373, 81463,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 114188, 0, 3,
                                                                       108868, 78628, 110128,
                                                                       57373, 58198, 82618,
                                                                       ncols, gamma, p, q);

                    compute_prim_msi_three_center_electron_repulsion_0(buffer, 115728, 0, 3,
                                                                       110128, 79573, 111388,
                                                                       58198, 59023, 83773,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117268, 3, 60673,
                                                                       60694, 84984, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117304, 3, 60694,
                                                                       60715, 85012, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117340, 3, 60715,
                                                                       60736, 85040, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117376, 3, 60736,
                                                                       60757, 85068, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117412, 3, 60757,
                                                                       60778, 85096, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117448, 3, 60778,
                                                                       60799, 85124, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117484, 3, 60799,
                                                                       60820, 85152, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117520, 3, 60820,
                                                                       60841, 85180, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117556, 3, 60841,
                                                                       60862, 85208, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 117592, 3, 60862,
                                                                       60883, 85236, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 117628, 0, 3,
                                                                       117268, 84984, 117304,
                                                                       60925, 60988, 85432,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 117736, 0, 3,
                                                                       117304, 85012, 117340,
                                                                       60988, 61051, 85516,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 117844, 0, 3,
                                                                       117340, 85040, 117376,
                                                                       61051, 61114, 85600,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 117952, 0, 3,
                                                                       117376, 85068, 117412,
                                                                       61114, 61177, 85684,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 118060, 0, 3,
                                                                       117412, 85096, 117448,
                                                                       61177, 61240, 85768,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 118168, 0, 3,
                                                                       117448, 85124, 117484,
                                                                       61240, 61303, 85852,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 118276, 0, 3,
                                                                       117484, 85152, 117520,
                                                                       61303, 61366, 85936,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 118384, 0, 3,
                                                                       117520, 85180, 117556,
                                                                       61366, 61429, 86020,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 118492, 0, 3,
                                                                       117556, 85208, 117592,
                                                                       61429, 61492, 86104,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 118600, 0, 3,
                                                                       117628, 85432, 117736,
                                                                       61618, 61744, 86524,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 118816, 0, 3,
                                                                       117736, 85516, 117844,
                                                                       61744, 61870, 86692,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 119032, 0, 3,
                                                                       117844, 85600, 117952,
                                                                       61870, 61996, 86860,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 119248, 0, 3,
                                                                       117952, 85684, 118060,
                                                                       61996, 62122, 87028,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 119464, 0, 3,
                                                                       118060, 85768, 118168,
                                                                       62122, 62248, 87196,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 119680, 0, 3,
                                                                       118168, 85852, 118276,
                                                                       62248, 62374, 87364,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 119896, 0, 3,
                                                                       118276, 85936, 118384,
                                                                       62374, 62500, 87532,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 120112, 0, 3,
                                                                       118384, 86020, 118492,
                                                                       62500, 62626, 87700,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 120328, 0, 3,
                                                                       118600, 86524, 118816,
                                                                       62878, 63088, 88428,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 120688, 0, 3,
                                                                       118816, 86692, 119032,
                                                                       63088, 63298, 88708,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 121048, 0, 3,
                                                                       119032, 86860, 119248,
                                                                       63298, 63508, 88988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 121408, 0, 3,
                                                                       119248, 87028, 119464,
                                                                       63508, 63718, 89268,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 121768, 0, 3,
                                                                       119464, 87196, 119680,
                                                                       63718, 63928, 89548,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 122128, 0, 3,
                                                                       119680, 87364, 119896,
                                                                       63928, 64138, 89828,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 122488, 0, 3,
                                                                       119896, 87532, 120112,
                                                                       64138, 64348, 90108,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 122848, 0, 3,
                                                                       120328, 88428, 120688,
                                                                       64768, 65083, 91228,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 123388, 0, 3,
                                                                       120688, 88708, 121048,
                                                                       65083, 65398, 91648,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 123928, 0, 3,
                                                                       121048, 88988, 121408,
                                                                       65398, 65713, 92068,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 124468, 0, 3,
                                                                       121408, 89268, 121768,
                                                                       65713, 66028, 92488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 125008, 0, 3,
                                                                       121768, 89548, 122128,
                                                                       66028, 66343, 92908,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 125548, 0, 3,
                                                                       122128, 89828, 122488,
                                                                       66343, 66658, 93328,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 126088, 0, 3,
                                                                       122848, 91228, 123388,
                                                                       67288, 67729, 94924,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 126844, 0, 3,
                                                                       123388, 91648, 123928,
                                                                       67729, 68170, 95512,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 127600, 0, 3,
                                                                       123928, 92068, 124468,
                                                                       68170, 68611, 96100,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 128356, 0, 3,
                                                                       124468, 92488, 125008,
                                                                       68611, 69052, 96688,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 129112, 0, 3,
                                                                       125008, 92908, 125548,
                                                                       69052, 69493, 97276,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 129868, 0, 3,
                                                                       126088, 94924, 126844,
                                                                       70375, 70963, 99432,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 130876, 0, 3,
                                                                       126844, 95512, 127600,
                                                                       70963, 71551, 100216,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 131884, 0, 3,
                                                                       127600, 96100, 128356,
                                                                       71551, 72139, 101000,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 132892, 0, 3,
                                                                       128356, 96688, 129112,
                                                                       72139, 72727, 101784,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 133900, 0, 3,
                                                                       129868, 99432, 130876,
                                                                       73903, 74659, 104584,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 135196, 0, 3,
                                                                       130876, 100216, 131884,
                                                                       74659, 75415, 105592,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 136492, 0, 3,
                                                                       131884, 101000, 132892,
                                                                       75415, 76171, 106600,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 137788, 0, 3,
                                                                       133900, 104584, 135196,
                                                                       77683, 78628, 110128,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 139408, 0, 3,
                                                                       135196, 105592, 136492,
                                                                       78628, 79573, 111388,
                                                                       ncols, gamma, p, q);

                    compute_prim_msk_three_center_electron_repulsion_0(buffer, 141028, 0, 3,
                                                                       137788, 110128, 139408,
                                                                       81463, 82618, 115728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143008, 3, 84928,
                                                                       84956, 117268, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143053, 3, 84956,
                                                                       84984, 117304, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143098, 3, 84984,
                                                                       85012, 117340, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143143, 3, 85012,
                                                                       85040, 117376, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143188, 3, 85040,
                                                                       85068, 117412, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143233, 3, 85068,
                                                                       85096, 117448, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143278, 3, 85096,
                                                                       85124, 117484, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143323, 3, 85124,
                                                                       85152, 117520, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143368, 3, 85152,
                                                                       85180, 117556, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 143413, 3, 85180,
                                                                       85208, 117592, ncols,
                                                                       gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 143458, 0, 3,
                                                                       143008, 117268, 143053,
                                                                       85264, 85348, 117628,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 143593, 0, 3,
                                                                       143053, 117304, 143098,
                                                                       85348, 85432, 117736,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 143728, 0, 3,
                                                                       143098, 117340, 143143,
                                                                       85432, 85516, 117844,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 143863, 0, 3,
                                                                       143143, 117376, 143188,
                                                                       85516, 85600, 117952,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 143998, 0, 3,
                                                                       143188, 117412, 143233,
                                                                       85600, 85684, 118060,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 144133, 0, 3,
                                                                       143233, 117448, 143278,
                                                                       85684, 85768, 118168,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 144268, 0, 3,
                                                                       143278, 117484, 143323,
                                                                       85768, 85852, 118276,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 144403, 0, 3,
                                                                       143323, 117520, 143368,
                                                                       85852, 85936, 118384,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 144538, 0, 3,
                                                                       143368, 117556, 143413,
                                                                       85936, 86020, 118492,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 144673, 0, 3,
                                                                       143458, 117628, 143593,
                                                                       86188, 86356, 118600,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 144943, 0, 3,
                                                                       143593, 117736, 143728,
                                                                       86356, 86524, 118816,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 145213, 0, 3,
                                                                       143728, 117844, 143863,
                                                                       86524, 86692, 119032,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 145483, 0, 3,
                                                                       143863, 117952, 143998,
                                                                       86692, 86860, 119248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 145753, 0, 3,
                                                                       143998, 118060, 144133,
                                                                       86860, 87028, 119464,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 146023, 0, 3,
                                                                       144133, 118168, 144268,
                                                                       87028, 87196, 119680,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 146293, 0, 3,
                                                                       144268, 118276, 144403,
                                                                       87196, 87364, 119896,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 146563, 0, 3,
                                                                       144403, 118384, 144538,
                                                                       87364, 87532, 120112,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 146833, 0, 3,
                                                                       144673, 118600, 144943,
                                                                       87868, 88148, 120328,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 147283, 0, 3,
                                                                       144943, 118816, 145213,
                                                                       88148, 88428, 120688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 147733, 0, 3,
                                                                       145213, 119032, 145483,
                                                                       88428, 88708, 121048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 148183, 0, 3,
                                                                       145483, 119248, 145753,
                                                                       88708, 88988, 121408,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 148633, 0, 3,
                                                                       145753, 119464, 146023,
                                                                       88988, 89268, 121768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 149083, 0, 3,
                                                                       146023, 119680, 146293,
                                                                       89268, 89548, 122128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 149533, 0, 3,
                                                                       146293, 119896, 146563,
                                                                       89548, 89828, 122488,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 149983, 0, 3,
                                                                       146833, 120328, 147283,
                                                                       90388, 90808, 122848,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 150658, 0, 3,
                                                                       147283, 120688, 147733,
                                                                       90808, 91228, 123388,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 151333, 0, 3,
                                                                       147733, 121048, 148183,
                                                                       91228, 91648, 123928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 152008, 0, 3,
                                                                       148183, 121408, 148633,
                                                                       91648, 92068, 124468,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 152683, 0, 3,
                                                                       148633, 121768, 149083,
                                                                       92068, 92488, 125008,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 153358, 0, 3,
                                                                       149083, 122128, 149533,
                                                                       92488, 92908, 125548,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 154033, 0, 3,
                                                                       149983, 122848, 150658,
                                                                       93748, 94336, 126088,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 154978, 0, 3,
                                                                       150658, 123388, 151333,
                                                                       94336, 94924, 126844,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 155923, 0, 3,
                                                                       151333, 123928, 152008,
                                                                       94924, 95512, 127600,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 156868, 0, 3,
                                                                       152008, 124468, 152683,
                                                                       95512, 96100, 128356,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 157813, 0, 3,
                                                                       152683, 125008, 153358,
                                                                       96100, 96688, 129112,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 158758, 0, 3,
                                                                       154033, 126088, 154978,
                                                                       97864, 98648, 129868,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 160018, 0, 3,
                                                                       154978, 126844, 155923,
                                                                       98648, 99432, 130876,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 161278, 0, 3,
                                                                       155923, 127600, 156868,
                                                                       99432, 100216, 131884,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 162538, 0, 3,
                                                                       156868, 128356, 157813,
                                                                       100216, 101000, 132892,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 163798, 0, 3,
                                                                       158758, 129868, 160018,
                                                                       102568, 103576, 133900,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 165418, 0, 3,
                                                                       160018, 130876, 161278,
                                                                       103576, 104584, 135196,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksl_three_center_electron_repulsion_0(buffer, 167038, 0, 3,
                                                                       161278, 131884, 162538,
                                                                       104584, 105592, 136492,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 168658, 0, 3,
                                                                       163798, 133900, 165418,
                                                                       107608, 108868, 137788,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsl_three_center_electron_repulsion_0(buffer, 170683, 0, 3,
                                                                       165418, 135196, 167038,
                                                                       108868, 110128, 139408,
                                                                       ncols, gamma, p, q);

                    compute_prim_msl_three_center_electron_repulsion_0(buffer, 172708, 0, 3,
                                                                       168658, 137788, 170683,
                                                                       112648, 114188, 141028,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 175183, 154033, 945, ncols);

                    simdfunc::contract_primitives(buffer, 176485, 158758, 1260, ncols);

                    simdfunc::contract_primitives(buffer, 178221, 163798, 1620, ncols);

                    simdfunc::contract_primitives(buffer, 180453, 168658, 2025, ncols);

                    simdfunc::contract_primitives(buffer, 183243, 172708, 2475, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 176128, 175183, 21, 1, nmax);

        simdtrf::transform_l_inner(buffer, 177745, 176485, 28, 1, nmax);

        simdtrf::transform_l_inner(buffer, 179841, 178221, 36, 1, nmax);

        simdtrf::transform_l_inner(buffer, 182478, 180453, 45, 1, nmax);

        simdtrf::transform_l_inner(buffer, 185718, 183243, 55, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 186653, 176128, 177745, 17,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 187724, 177745, 179841, 17,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 189152, 179841, 182478, 17,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 190988, 182478, 185718, 17,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 193283, 186653, 187724, 17,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 195425, 187724, 189152, 17,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 198281, 189152, 190988, 17,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 201953, 193283, 195425, 17,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 205523, 195425, 198281, 17,
                                             nmax);

        simdtrf::compute_hrr_hg_out_of_first(buffer, coordinates, 210283, 201953, 205523, 17,
                                             nmax);

        simdtrf::transform_g_inner(buffer, 215638, 210283, 21, 17, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 215638, 153, nmax);
    }

    for (size_t m = 0; m < 1683; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
