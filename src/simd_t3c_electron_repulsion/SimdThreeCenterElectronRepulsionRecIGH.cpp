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


#include "SimdThreeCenterElectronRepulsionRecIGH.hpp"

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
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_igh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_igh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 104189, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 104189, 69791, 6634, dimensions);

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
                                                        5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
                                                        ncols, fj, mu, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 22, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

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

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 7, 8,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 8, 9,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 76, 0, 3, 9, 10,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 88, 0, 3, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 94, 0, 3, 12, 13,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 100, 0, 3, 13, 14,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 106, 0, 3, 14, 15,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 15, 16,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 16, 17,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 124, 0, 3, 17, 18,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 130, 0, 3, 18, 19,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 136, 0, 3, 19, 20,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 22, 25,
                                                                       64, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 25, 28,
                                                                       70, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 162, 0, 3, 28, 31,
                                                                       76, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 31, 34,
                                                                       82, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 34, 37,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 37, 40,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 40, 43,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 212, 0, 3, 43, 46,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 222, 0, 3, 46, 49,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 49, 52,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 242, 0, 3, 52, 55,
                                                                       124, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 252, 0, 3, 55, 58,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 262, 0, 3, 64, 70,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 277, 0, 3, 70, 76,
                                                                       152, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 292, 0, 3, 76, 82,
                                                                       162, 172, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 307, 0, 3, 82, 88,
                                                                       172, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 322, 0, 3, 88, 94,
                                                                       182, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 337, 0, 3, 94,
                                                                       100, 192, 202, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 352, 0, 3, 100,
                                                                       106, 202, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 367, 0, 3, 106,
                                                                       112, 212, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 382, 0, 3, 112,
                                                                       118, 222, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 397, 0, 3, 118,
                                                                       124, 232, 242, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 412, 0, 3, 124,
                                                                       130, 242, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 427, 0, 3, 142,
                                                                       152, 262, 277, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 152,
                                                                       162, 277, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 469, 0, 3, 162,
                                                                       172, 292, 307, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 490, 0, 3, 172,
                                                                       182, 307, 322, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 511, 0, 3, 182,
                                                                       192, 322, 337, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 532, 0, 3, 192,
                                                                       202, 337, 352, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 553, 0, 3, 202,
                                                                       212, 352, 367, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 574, 0, 3, 212,
                                                                       222, 367, 382, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 595, 0, 3, 222,
                                                                       232, 382, 397, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 616, 0, 3, 232,
                                                                       242, 397, 412, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 637, 0, 3, 262,
                                                                       277, 427, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 665, 0, 3, 277,
                                                                       292, 448, 469, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 693, 0, 3, 292,
                                                                       307, 469, 490, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 721, 0, 3, 307,
                                                                       322, 490, 511, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 749, 0, 3, 322,
                                                                       337, 511, 532, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 777, 0, 3, 337,
                                                                       352, 532, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 805, 0, 3, 352,
                                                                       367, 553, 574, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 833, 0, 3, 367,
                                                                       382, 574, 595, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 861, 0, 3, 382,
                                                                       397, 595, 616, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 889, 0, 3, 427,
                                                                       448, 637, 665, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 925, 0, 3, 448,
                                                                       469, 665, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 961, 0, 3, 469,
                                                                       490, 693, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 997, 0, 3, 490,
                                                                       511, 721, 749, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1033, 0, 3, 511,
                                                                       532, 749, 777, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1069, 0, 3, 532,
                                                                       553, 777, 805, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1105, 0, 3, 553,
                                                                       574, 805, 833, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1141, 0, 3, 574,
                                                                       595, 833, 861, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1177, 0, 3, 637,
                                                                       665, 889, 925, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1222, 0, 3, 665,
                                                                       693, 925, 961, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1267, 0, 3, 693,
                                                                       721, 961, 997, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1312, 0, 3, 721,
                                                                       749, 997, 1033, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1357, 0, 3, 749,
                                                                       777, 1033, 1069, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1402, 0, 3, 777,
                                                                       805, 1069, 1105, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1447, 0, 3, 805,
                                                                       833, 1105, 1141, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 889,
                                                                       925, 1177, 1222, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1547, 0, 3, 925,
                                                                       961, 1222, 1267, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1602, 0, 3, 961,
                                                                       997, 1267, 1312, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1657, 0, 3, 997,
                                                                       1033, 1312, 1357, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1712, 0, 3, 1033,
                                                                       1069, 1357, 1402, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 1767, 0, 3, 1069,
                                                                       1105, 1402, 1447, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1822, 0, 3, 1177,
                                                                       1222, 1492, 1547, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1888, 0, 3, 1222,
                                                                       1267, 1547, 1602, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 1954, 0, 3, 1267,
                                                                       1312, 1602, 1657, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2020, 0, 3, 1312,
                                                                       1357, 1657, 1712, ncols,
                                                                       gamma, p, q);

                    compute_prim_nss_three_center_electron_repulsion_0(buffer, 2086, 0, 3, 1357,
                                                                       1402, 1712, 1767, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2152, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2155, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2158, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2161, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2164, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2167, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2170, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2173, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2176, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2179, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2182, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2185, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2188, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2191, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2194, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2197, 3, 7, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2206, 3, 8, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2215, 3, 9, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2224, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2233, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2242, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2251, 3, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2260, 3, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2269, 3, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2278, 3, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2287, 3, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2296, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2305, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2314, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2323, 3, 22, 64,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2341, 3, 25, 70,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2359, 3, 28, 76,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2377, 3, 31, 82,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2395, 3, 34, 88,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2413, 3, 37, 94,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2431, 3, 40, 100,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2449, 3, 43, 106,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2467, 3, 46, 112,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2485, 3, 49, 118,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2503, 3, 52, 124,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2521, 3, 55, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2539, 3, 58, 136,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2557, 3, 64, 142,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2587, 3, 70, 152,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2617, 3, 76, 162,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2647, 3, 82, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2677, 3, 88, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2707, 3, 94, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2737, 3, 100, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2767, 3, 106, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2797, 3, 112, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2827, 3, 118, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2857, 3, 124, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2887, 3, 130, 252,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2917, 3, 142, 262,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2962, 3, 152, 277,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3007, 3, 162, 292,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3052, 3, 172, 307,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3097, 3, 182, 322,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3142, 3, 192, 337,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3187, 3, 202, 352,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3232, 3, 212, 367,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3277, 3, 222, 382,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3322, 3, 232, 397,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3367, 3, 242, 412,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3412, 3, 262, 427,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3475, 3, 277, 448,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3538, 3, 292, 469,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3601, 3, 307, 490,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3664, 3, 322, 511,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3727, 3, 337, 532,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3790, 3, 352, 553,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3853, 3, 367, 574,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3916, 3, 382, 595,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3979, 3, 397, 616,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4042, 3, 427, 637,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4126, 3, 448, 665,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4210, 3, 469, 693,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4294, 3, 490, 721,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4378, 3, 511, 749,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4462, 3, 532, 777,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4546, 3, 553, 805,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4630, 3, 574, 833,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4714, 3, 595, 861,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4798, 3, 637, 889,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4906, 3, 665, 925,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5014, 3, 693, 961,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5122, 3, 721, 997,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5230, 3, 749,
                                                                       1033, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5338, 3, 777,
                                                                       1069, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5446, 3, 805,
                                                                       1105, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5554, 3, 833,
                                                                       1141, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5662, 3, 889,
                                                                       1177, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5797, 3, 925,
                                                                       1222, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5932, 3, 961,
                                                                       1267, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6067, 3, 997,
                                                                       1312, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6202, 3, 1033,
                                                                       1357, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6337, 3, 1069,
                                                                       1402, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6472, 3, 1105,
                                                                       1447, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6607, 3, 1177,
                                                                       1492, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6772, 3, 1222,
                                                                       1547, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 6937, 3, 1267,
                                                                       1602, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7102, 3, 1312,
                                                                       1657, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7267, 3, 1357,
                                                                       1712, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7432, 3, 1402,
                                                                       1767, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 7597, 3, 1492,
                                                                       1822, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 7795, 3, 1547,
                                                                       1888, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 7993, 3, 1602,
                                                                       1954, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8191, 3, 1657,
                                                                       2020, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8389, 3, 1712,
                                                                       2086, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8587, 3, 7, 8,
                                                                       2158, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8593, 3, 8, 9,
                                                                       2161, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8599, 3, 9, 10,
                                                                       2164, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8605, 3, 10, 11,
                                                                       2167, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8611, 3, 11, 12,
                                                                       2170, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8617, 3, 12, 13,
                                                                       2173, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8623, 3, 13, 14,
                                                                       2176, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8629, 3, 14, 15,
                                                                       2179, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8635, 3, 15, 16,
                                                                       2182, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8641, 3, 16, 17,
                                                                       2185, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8647, 3, 17, 18,
                                                                       2188, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8653, 3, 18, 19,
                                                                       2191, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 8659, 3, 19, 20,
                                                                       2194, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8665, 0, 3, 8587,
                                                                       2158, 8593, 2215, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8683, 0, 3, 8593,
                                                                       2161, 8599, 2224, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8701, 0, 3, 8599,
                                                                       2164, 8605, 2233, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8719, 0, 3, 8605,
                                                                       2167, 8611, 2242, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8737, 0, 3, 8611,
                                                                       2170, 8617, 2251, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8755, 0, 3, 8617,
                                                                       2173, 8623, 2260, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8773, 0, 3, 8623,
                                                                       2176, 8629, 2269, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8791, 0, 3, 8629,
                                                                       2179, 8635, 2278, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8809, 0, 3, 8635,
                                                                       2182, 8641, 2287, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8827, 0, 3, 8641,
                                                                       2185, 8647, 2296, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8845, 0, 3, 8647,
                                                                       2188, 8653, 2305, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8863, 0, 3, 8653,
                                                                       2191, 8659, 2314, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8881, 0, 3, 8665,
                                                                       2215, 8683, 64, 70, 2359,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8917, 0, 3, 8683,
                                                                       2224, 8701, 70, 76, 2377,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8953, 0, 3, 8701,
                                                                       2233, 8719, 76, 82, 2395,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8989, 0, 3, 8719,
                                                                       2242, 8737, 82, 88, 2413,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9025, 0, 3, 8737,
                                                                       2251, 8755, 88, 94, 2431,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9061, 0, 3, 8755,
                                                                       2260, 8773, 94, 100, 2449,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9097, 0, 3, 8773,
                                                                       2269, 8791, 100, 106,
                                                                       2467, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9133, 0, 3, 8791,
                                                                       2278, 8809, 106, 112,
                                                                       2485, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9169, 0, 3, 8809,
                                                                       2287, 8827, 112, 118,
                                                                       2503, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9205, 0, 3, 8827,
                                                                       2296, 8845, 118, 124,
                                                                       2521, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9241, 0, 3, 8845,
                                                                       2305, 8863, 124, 130,
                                                                       2539, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9277, 0, 3, 8881,
                                                                       2359, 8917, 142, 152,
                                                                       2617, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9337, 0, 3, 8917,
                                                                       2377, 8953, 152, 162,
                                                                       2647, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9397, 0, 3, 8953,
                                                                       2395, 8989, 162, 172,
                                                                       2677, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9457, 0, 3, 8989,
                                                                       2413, 9025, 172, 182,
                                                                       2707, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9517, 0, 3, 9025,
                                                                       2431, 9061, 182, 192,
                                                                       2737, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9577, 0, 3, 9061,
                                                                       2449, 9097, 192, 202,
                                                                       2767, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9637, 0, 3, 9097,
                                                                       2467, 9133, 202, 212,
                                                                       2797, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9697, 0, 3, 9133,
                                                                       2485, 9169, 212, 222,
                                                                       2827, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9757, 0, 3, 9169,
                                                                       2503, 9205, 222, 232,
                                                                       2857, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9817, 0, 3, 9205,
                                                                       2521, 9241, 232, 242,
                                                                       2887, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9877, 0, 3, 9277,
                                                                       2617, 9337, 262, 277,
                                                                       3007, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9967, 0, 3, 9337,
                                                                       2647, 9397, 277, 292,
                                                                       3052, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10057, 0, 3, 9397,
                                                                       2677, 9457, 292, 307,
                                                                       3097, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10147, 0, 3, 9457,
                                                                       2707, 9517, 307, 322,
                                                                       3142, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10237, 0, 3, 9517,
                                                                       2737, 9577, 322, 337,
                                                                       3187, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10327, 0, 3, 9577,
                                                                       2767, 9637, 337, 352,
                                                                       3232, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10417, 0, 3, 9637,
                                                                       2797, 9697, 352, 367,
                                                                       3277, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10507, 0, 3, 9697,
                                                                       2827, 9757, 367, 382,
                                                                       3322, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10597, 0, 3, 9757,
                                                                       2857, 9817, 382, 397,
                                                                       3367, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10687, 0, 3, 9877,
                                                                       3007, 9967, 427, 448,
                                                                       3538, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10813, 0, 3, 9967,
                                                                       3052, 10057, 448, 469,
                                                                       3601, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10939, 0, 3,
                                                                       10057, 3097, 10147, 469,
                                                                       490, 3664, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11065, 0, 3,
                                                                       10147, 3142, 10237, 490,
                                                                       511, 3727, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11191, 0, 3,
                                                                       10237, 3187, 10327, 511,
                                                                       532, 3790, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11317, 0, 3,
                                                                       10327, 3232, 10417, 532,
                                                                       553, 3853, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11443, 0, 3,
                                                                       10417, 3277, 10507, 553,
                                                                       574, 3916, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11569, 0, 3,
                                                                       10507, 3322, 10597, 574,
                                                                       595, 3979, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11695, 0, 3,
                                                                       10687, 3538, 10813, 637,
                                                                       665, 4210, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 11863, 0, 3,
                                                                       10813, 3601, 10939, 665,
                                                                       693, 4294, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12031, 0, 3,
                                                                       10939, 3664, 11065, 693,
                                                                       721, 4378, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12199, 0, 3,
                                                                       11065, 3727, 11191, 721,
                                                                       749, 4462, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12367, 0, 3,
                                                                       11191, 3790, 11317, 749,
                                                                       777, 4546, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12535, 0, 3,
                                                                       11317, 3853, 11443, 777,
                                                                       805, 4630, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12703, 0, 3,
                                                                       11443, 3916, 11569, 805,
                                                                       833, 4714, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 12871, 0, 3,
                                                                       11695, 4210, 11863, 889,
                                                                       925, 5014, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13087, 0, 3,
                                                                       11863, 4294, 12031, 925,
                                                                       961, 5122, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13303, 0, 3,
                                                                       12031, 4378, 12199, 961,
                                                                       997, 5230, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13519, 0, 3,
                                                                       12199, 4462, 12367, 997,
                                                                       1033, 5338, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13735, 0, 3,
                                                                       12367, 4546, 12535, 1033,
                                                                       1069, 5446, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 13951, 0, 3,
                                                                       12535, 4630, 12703, 1069,
                                                                       1105, 5554, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14167, 0, 3,
                                                                       12871, 5014, 13087, 1177,
                                                                       1222, 5932, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14437, 0, 3,
                                                                       13087, 5122, 13303, 1222,
                                                                       1267, 6067, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14707, 0, 3,
                                                                       13303, 5230, 13519, 1267,
                                                                       1312, 6202, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 14977, 0, 3,
                                                                       13519, 5338, 13735, 1312,
                                                                       1357, 6337, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 15247, 0, 3,
                                                                       13735, 5446, 13951, 1357,
                                                                       1402, 6472, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 15517, 0, 3,
                                                                       14167, 5932, 14437, 1492,
                                                                       1547, 6937, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 15847, 0, 3,
                                                                       14437, 6067, 14707, 1547,
                                                                       1602, 7102, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 16177, 0, 3,
                                                                       14707, 6202, 14977, 1602,
                                                                       1657, 7267, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 16507, 0, 3,
                                                                       14977, 6337, 15247, 1657,
                                                                       1712, 7432, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 16837, 0, 3,
                                                                       15517, 6937, 15847, 1822,
                                                                       1888, 7993, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 17233, 0, 3,
                                                                       15847, 7102, 16177, 1888,
                                                                       1954, 8191, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 17629, 0, 3,
                                                                       16177, 7267, 16507, 1954,
                                                                       2020, 8389, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18025, 3, 2152,
                                                                       2155, 8587, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18035, 3, 2155,
                                                                       2158, 8593, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18045, 3, 2158,
                                                                       2161, 8599, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18055, 3, 2161,
                                                                       2164, 8605, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18065, 3, 2164,
                                                                       2167, 8611, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18075, 3, 2167,
                                                                       2170, 8617, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18085, 3, 2170,
                                                                       2173, 8623, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18095, 3, 2173,
                                                                       2176, 8629, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18105, 3, 2176,
                                                                       2179, 8635, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18115, 3, 2179,
                                                                       2182, 8641, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18125, 3, 2182,
                                                                       2185, 8647, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18135, 3, 2185,
                                                                       2188, 8653, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18145, 3, 2188,
                                                                       2191, 8659, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18155, 0, 3,
                                                                       18025, 8587, 18035, 2197,
                                                                       2206, 8665, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18185, 0, 3,
                                                                       18035, 8593, 18045, 2206,
                                                                       2215, 8683, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18215, 0, 3,
                                                                       18045, 8599, 18055, 2215,
                                                                       2224, 8701, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18245, 0, 3,
                                                                       18055, 8605, 18065, 2224,
                                                                       2233, 8719, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18275, 0, 3,
                                                                       18065, 8611, 18075, 2233,
                                                                       2242, 8737, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18305, 0, 3,
                                                                       18075, 8617, 18085, 2242,
                                                                       2251, 8755, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18335, 0, 3,
                                                                       18085, 8623, 18095, 2251,
                                                                       2260, 8773, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18365, 0, 3,
                                                                       18095, 8629, 18105, 2260,
                                                                       2269, 8791, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18395, 0, 3,
                                                                       18105, 8635, 18115, 2269,
                                                                       2278, 8809, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18425, 0, 3,
                                                                       18115, 8641, 18125, 2278,
                                                                       2287, 8827, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18455, 0, 3,
                                                                       18125, 8647, 18135, 2287,
                                                                       2296, 8845, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18485, 0, 3,
                                                                       18135, 8653, 18145, 2296,
                                                                       2305, 8863, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18515, 0, 3,
                                                                       18155, 8665, 18185, 2323,
                                                                       2341, 8881, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18575, 0, 3,
                                                                       18185, 8683, 18215, 2341,
                                                                       2359, 8917, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18635, 0, 3,
                                                                       18215, 8701, 18245, 2359,
                                                                       2377, 8953, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18695, 0, 3,
                                                                       18245, 8719, 18275, 2377,
                                                                       2395, 8989, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18755, 0, 3,
                                                                       18275, 8737, 18305, 2395,
                                                                       2413, 9025, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18815, 0, 3,
                                                                       18305, 8755, 18335, 2413,
                                                                       2431, 9061, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18875, 0, 3,
                                                                       18335, 8773, 18365, 2431,
                                                                       2449, 9097, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18935, 0, 3,
                                                                       18365, 8791, 18395, 2449,
                                                                       2467, 9133, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18995, 0, 3,
                                                                       18395, 8809, 18425, 2467,
                                                                       2485, 9169, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19055, 0, 3,
                                                                       18425, 8827, 18455, 2485,
                                                                       2503, 9205, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 19115, 0, 3,
                                                                       18455, 8845, 18485, 2503,
                                                                       2521, 9241, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19175, 0, 3,
                                                                       18515, 8881, 18575, 2557,
                                                                       2587, 9277, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19275, 0, 3,
                                                                       18575, 8917, 18635, 2587,
                                                                       2617, 9337, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19375, 0, 3,
                                                                       18635, 8953, 18695, 2617,
                                                                       2647, 9397, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19475, 0, 3,
                                                                       18695, 8989, 18755, 2647,
                                                                       2677, 9457, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19575, 0, 3,
                                                                       18755, 9025, 18815, 2677,
                                                                       2707, 9517, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19675, 0, 3,
                                                                       18815, 9061, 18875, 2707,
                                                                       2737, 9577, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19775, 0, 3,
                                                                       18875, 9097, 18935, 2737,
                                                                       2767, 9637, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19875, 0, 3,
                                                                       18935, 9133, 18995, 2767,
                                                                       2797, 9697, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19975, 0, 3,
                                                                       18995, 9169, 19055, 2797,
                                                                       2827, 9757, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 20075, 0, 3,
                                                                       19055, 9205, 19115, 2827,
                                                                       2857, 9817, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20175, 0, 3,
                                                                       19175, 9277, 19275, 2917,
                                                                       2962, 9877, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20325, 0, 3,
                                                                       19275, 9337, 19375, 2962,
                                                                       3007, 9967, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20475, 0, 3,
                                                                       19375, 9397, 19475, 3007,
                                                                       3052, 10057, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20625, 0, 3,
                                                                       19475, 9457, 19575, 3052,
                                                                       3097, 10147, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20775, 0, 3,
                                                                       19575, 9517, 19675, 3097,
                                                                       3142, 10237, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20925, 0, 3,
                                                                       19675, 9577, 19775, 3142,
                                                                       3187, 10327, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21075, 0, 3,
                                                                       19775, 9637, 19875, 3187,
                                                                       3232, 10417, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21225, 0, 3,
                                                                       19875, 9697, 19975, 3232,
                                                                       3277, 10507, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21375, 0, 3,
                                                                       19975, 9757, 20075, 3277,
                                                                       3322, 10597, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 21525, 0, 3,
                                                                       20175, 9877, 20325, 3412,
                                                                       3475, 10687, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 21735, 0, 3,
                                                                       20325, 9967, 20475, 3475,
                                                                       3538, 10813, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 21945, 0, 3,
                                                                       20475, 10057, 20625, 3538,
                                                                       3601, 10939, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22155, 0, 3,
                                                                       20625, 10147, 20775, 3601,
                                                                       3664, 11065, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22365, 0, 3,
                                                                       20775, 10237, 20925, 3664,
                                                                       3727, 11191, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22575, 0, 3,
                                                                       20925, 10327, 21075, 3727,
                                                                       3790, 11317, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22785, 0, 3,
                                                                       21075, 10417, 21225, 3790,
                                                                       3853, 11443, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22995, 0, 3,
                                                                       21225, 10507, 21375, 3853,
                                                                       3916, 11569, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 23205, 0, 3,
                                                                       21525, 10687, 21735, 4042,
                                                                       4126, 11695, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 23485, 0, 3,
                                                                       21735, 10813, 21945, 4126,
                                                                       4210, 11863, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 23765, 0, 3,
                                                                       21945, 10939, 22155, 4210,
                                                                       4294, 12031, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24045, 0, 3,
                                                                       22155, 11065, 22365, 4294,
                                                                       4378, 12199, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24325, 0, 3,
                                                                       22365, 11191, 22575, 4378,
                                                                       4462, 12367, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24605, 0, 3,
                                                                       22575, 11317, 22785, 4462,
                                                                       4546, 12535, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24885, 0, 3,
                                                                       22785, 11443, 22995, 4546,
                                                                       4630, 12703, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 25165, 0, 3,
                                                                       23205, 11695, 23485, 4798,
                                                                       4906, 12871, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 25525, 0, 3,
                                                                       23485, 11863, 23765, 4906,
                                                                       5014, 13087, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 25885, 0, 3,
                                                                       23765, 12031, 24045, 5014,
                                                                       5122, 13303, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 26245, 0, 3,
                                                                       24045, 12199, 24325, 5122,
                                                                       5230, 13519, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 26605, 0, 3,
                                                                       24325, 12367, 24605, 5230,
                                                                       5338, 13735, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 26965, 0, 3,
                                                                       24605, 12535, 24885, 5338,
                                                                       5446, 13951, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 27325, 0, 3,
                                                                       25165, 12871, 25525, 5662,
                                                                       5797, 14167, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 27775, 0, 3,
                                                                       25525, 13087, 25885, 5797,
                                                                       5932, 14437, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 28225, 0, 3,
                                                                       25885, 13303, 26245, 5932,
                                                                       6067, 14707, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 28675, 0, 3,
                                                                       26245, 13519, 26605, 6067,
                                                                       6202, 14977, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 29125, 0, 3,
                                                                       26605, 13735, 26965, 6202,
                                                                       6337, 15247, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 29575, 0, 3,
                                                                       27325, 14167, 27775, 6607,
                                                                       6772, 15517, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 30125, 0, 3,
                                                                       27775, 14437, 28225, 6772,
                                                                       6937, 15847, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 30675, 0, 3,
                                                                       28225, 14707, 28675, 6937,
                                                                       7102, 16177, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 31225, 0, 3,
                                                                       28675, 14977, 29125, 7102,
                                                                       7267, 16507, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 31775, 0, 3,
                                                                       29575, 15517, 30125, 7597,
                                                                       7795, 16837, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 32435, 0, 3,
                                                                       30125, 15847, 30675, 7795,
                                                                       7993, 17233, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 33095, 0, 3,
                                                                       30675, 16177, 31225, 7993,
                                                                       8191, 17629, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33755, 3, 8587,
                                                                       8593, 18045, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33770, 3, 8593,
                                                                       8599, 18055, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33785, 3, 8599,
                                                                       8605, 18065, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33800, 3, 8605,
                                                                       8611, 18075, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33815, 3, 8611,
                                                                       8617, 18085, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33830, 3, 8617,
                                                                       8623, 18095, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33845, 3, 8623,
                                                                       8629, 18105, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33860, 3, 8629,
                                                                       8635, 18115, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33875, 3, 8635,
                                                                       8641, 18125, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33890, 3, 8641,
                                                                       8647, 18135, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 33905, 3, 8647,
                                                                       8653, 18145, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 33920, 0, 3,
                                                                       33755, 18045, 33770, 8665,
                                                                       8683, 18215, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 33965, 0, 3,
                                                                       33770, 18055, 33785, 8683,
                                                                       8701, 18245, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34010, 0, 3,
                                                                       33785, 18065, 33800, 8701,
                                                                       8719, 18275, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34055, 0, 3,
                                                                       33800, 18075, 33815, 8719,
                                                                       8737, 18305, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34100, 0, 3,
                                                                       33815, 18085, 33830, 8737,
                                                                       8755, 18335, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34145, 0, 3,
                                                                       33830, 18095, 33845, 8755,
                                                                       8773, 18365, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34190, 0, 3,
                                                                       33845, 18105, 33860, 8773,
                                                                       8791, 18395, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34235, 0, 3,
                                                                       33860, 18115, 33875, 8791,
                                                                       8809, 18425, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34280, 0, 3,
                                                                       33875, 18125, 33890, 8809,
                                                                       8827, 18455, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 34325, 0, 3,
                                                                       33890, 18135, 33905, 8827,
                                                                       8845, 18485, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34370, 0, 3,
                                                                       33920, 18215, 33965, 8881,
                                                                       8917, 18635, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34460, 0, 3,
                                                                       33965, 18245, 34010, 8917,
                                                                       8953, 18695, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34550, 0, 3,
                                                                       34010, 18275, 34055, 8953,
                                                                       8989, 18755, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34640, 0, 3,
                                                                       34055, 18305, 34100, 8989,
                                                                       9025, 18815, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34730, 0, 3,
                                                                       34100, 18335, 34145, 9025,
                                                                       9061, 18875, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34820, 0, 3,
                                                                       34145, 18365, 34190, 9061,
                                                                       9097, 18935, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34910, 0, 3,
                                                                       34190, 18395, 34235, 9097,
                                                                       9133, 18995, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35000, 0, 3,
                                                                       34235, 18425, 34280, 9133,
                                                                       9169, 19055, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 35090, 0, 3,
                                                                       34280, 18455, 34325, 9169,
                                                                       9205, 19115, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35180, 0, 3,
                                                                       34370, 18635, 34460, 9277,
                                                                       9337, 19375, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35330, 0, 3,
                                                                       34460, 18695, 34550, 9337,
                                                                       9397, 19475, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35480, 0, 3,
                                                                       34550, 18755, 34640, 9397,
                                                                       9457, 19575, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35630, 0, 3,
                                                                       34640, 18815, 34730, 9457,
                                                                       9517, 19675, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35780, 0, 3,
                                                                       34730, 18875, 34820, 9517,
                                                                       9577, 19775, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 35930, 0, 3,
                                                                       34820, 18935, 34910, 9577,
                                                                       9637, 19875, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 36080, 0, 3,
                                                                       34910, 18995, 35000, 9637,
                                                                       9697, 19975, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 36230, 0, 3,
                                                                       35000, 19055, 35090, 9697,
                                                                       9757, 20075, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 36380, 0, 3,
                                                                       35180, 19375, 35330, 9877,
                                                                       9967, 20475, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 36605, 0, 3,
                                                                       35330, 19475, 35480, 9967,
                                                                       10057, 20625, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 36830, 0, 3,
                                                                       35480, 19575, 35630,
                                                                       10057, 10147, 20775,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37055, 0, 3,
                                                                       35630, 19675, 35780,
                                                                       10147, 10237, 20925,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37280, 0, 3,
                                                                       35780, 19775, 35930,
                                                                       10237, 10327, 21075,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37505, 0, 3,
                                                                       35930, 19875, 36080,
                                                                       10327, 10417, 21225,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 37730, 0, 3,
                                                                       36080, 19975, 36230,
                                                                       10417, 10507, 21375,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 37955, 0, 3,
                                                                       36380, 20475, 36605,
                                                                       10687, 10813, 21945,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 38270, 0, 3,
                                                                       36605, 20625, 36830,
                                                                       10813, 10939, 22155,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 38585, 0, 3,
                                                                       36830, 20775, 37055,
                                                                       10939, 11065, 22365,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 38900, 0, 3,
                                                                       37055, 20925, 37280,
                                                                       11065, 11191, 22575,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 39215, 0, 3,
                                                                       37280, 21075, 37505,
                                                                       11191, 11317, 22785,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 39530, 0, 3,
                                                                       37505, 21225, 37730,
                                                                       11317, 11443, 22995,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 39845, 0, 3,
                                                                       37955, 21945, 38270,
                                                                       11695, 11863, 23765,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 40265, 0, 3,
                                                                       38270, 22155, 38585,
                                                                       11863, 12031, 24045,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 40685, 0, 3,
                                                                       38585, 22365, 38900,
                                                                       12031, 12199, 24325,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 41105, 0, 3,
                                                                       38900, 22575, 39215,
                                                                       12199, 12367, 24605,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 41525, 0, 3,
                                                                       39215, 22785, 39530,
                                                                       12367, 12535, 24885,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 41945, 0, 3,
                                                                       39845, 23765, 40265,
                                                                       12871, 13087, 25885,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 42485, 0, 3,
                                                                       40265, 24045, 40685,
                                                                       13087, 13303, 26245,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 43025, 0, 3,
                                                                       40685, 24325, 41105,
                                                                       13303, 13519, 26605,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 43565, 0, 3,
                                                                       41105, 24605, 41525,
                                                                       13519, 13735, 26965,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 44105, 0, 3,
                                                                       41945, 25885, 42485,
                                                                       14167, 14437, 28225,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 44780, 0, 3,
                                                                       42485, 26245, 43025,
                                                                       14437, 14707, 28675,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 45455, 0, 3,
                                                                       43025, 26605, 43565,
                                                                       14707, 14977, 29125,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 46130, 0, 3,
                                                                       44105, 28225, 44780,
                                                                       15517, 15847, 30675,
                                                                       ncols, gamma, p, q);

                    compute_prim_msg_three_center_electron_repulsion_0(buffer, 46955, 0, 3,
                                                                       44780, 28675, 45455,
                                                                       15847, 16177, 31225,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsg_three_center_electron_repulsion_0(buffer, 47780, 0, 3,
                                                                       46130, 30675, 46955,
                                                                       16837, 17233, 33095,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48770, 3, 18025,
                                                                       18035, 33755, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48791, 3, 18035,
                                                                       18045, 33770, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48812, 3, 18045,
                                                                       18055, 33785, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48833, 3, 18055,
                                                                       18065, 33800, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48854, 3, 18065,
                                                                       18075, 33815, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48875, 3, 18075,
                                                                       18085, 33830, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48896, 3, 18085,
                                                                       18095, 33845, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48917, 3, 18095,
                                                                       18105, 33860, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48938, 3, 18105,
                                                                       18115, 33875, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48959, 3, 18115,
                                                                       18125, 33890, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48980, 3, 18125,
                                                                       18135, 33905, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49001, 0, 3,
                                                                       48770, 33755, 48791,
                                                                       18155, 18185, 33920,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49064, 0, 3,
                                                                       48791, 33770, 48812,
                                                                       18185, 18215, 33965,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49127, 0, 3,
                                                                       48812, 33785, 48833,
                                                                       18215, 18245, 34010,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49190, 0, 3,
                                                                       48833, 33800, 48854,
                                                                       18245, 18275, 34055,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49253, 0, 3,
                                                                       48854, 33815, 48875,
                                                                       18275, 18305, 34100,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49316, 0, 3,
                                                                       48875, 33830, 48896,
                                                                       18305, 18335, 34145,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49379, 0, 3,
                                                                       48896, 33845, 48917,
                                                                       18335, 18365, 34190,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49442, 0, 3,
                                                                       48917, 33860, 48938,
                                                                       18365, 18395, 34235,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49505, 0, 3,
                                                                       48938, 33875, 48959,
                                                                       18395, 18425, 34280,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49568, 0, 3,
                                                                       48959, 33890, 48980,
                                                                       18425, 18455, 34325,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 49631, 0, 3,
                                                                       49001, 33920, 49064,
                                                                       18515, 18575, 34370,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 49757, 0, 3,
                                                                       49064, 33965, 49127,
                                                                       18575, 18635, 34460,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 49883, 0, 3,
                                                                       49127, 34010, 49190,
                                                                       18635, 18695, 34550,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 50009, 0, 3,
                                                                       49190, 34055, 49253,
                                                                       18695, 18755, 34640,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 50135, 0, 3,
                                                                       49253, 34100, 49316,
                                                                       18755, 18815, 34730,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 50261, 0, 3,
                                                                       49316, 34145, 49379,
                                                                       18815, 18875, 34820,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 50387, 0, 3,
                                                                       49379, 34190, 49442,
                                                                       18875, 18935, 34910,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 50513, 0, 3,
                                                                       49442, 34235, 49505,
                                                                       18935, 18995, 35000,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 50639, 0, 3,
                                                                       49505, 34280, 49568,
                                                                       18995, 19055, 35090,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 50765, 0, 3,
                                                                       49631, 34370, 49757,
                                                                       19175, 19275, 35180,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 50975, 0, 3,
                                                                       49757, 34460, 49883,
                                                                       19275, 19375, 35330,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 51185, 0, 3,
                                                                       49883, 34550, 50009,
                                                                       19375, 19475, 35480,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 51395, 0, 3,
                                                                       50009, 34640, 50135,
                                                                       19475, 19575, 35630,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 51605, 0, 3,
                                                                       50135, 34730, 50261,
                                                                       19575, 19675, 35780,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 51815, 0, 3,
                                                                       50261, 34820, 50387,
                                                                       19675, 19775, 35930,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 52025, 0, 3,
                                                                       50387, 34910, 50513,
                                                                       19775, 19875, 36080,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 52235, 0, 3,
                                                                       50513, 35000, 50639,
                                                                       19875, 19975, 36230,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 52445, 0, 3,
                                                                       50765, 35180, 50975,
                                                                       20175, 20325, 36380,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 52760, 0, 3,
                                                                       50975, 35330, 51185,
                                                                       20325, 20475, 36605,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 53075, 0, 3,
                                                                       51185, 35480, 51395,
                                                                       20475, 20625, 36830,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 53390, 0, 3,
                                                                       51395, 35630, 51605,
                                                                       20625, 20775, 37055,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 53705, 0, 3,
                                                                       51605, 35780, 51815,
                                                                       20775, 20925, 37280,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 54020, 0, 3,
                                                                       51815, 35930, 52025,
                                                                       20925, 21075, 37505,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 54335, 0, 3,
                                                                       52025, 36080, 52235,
                                                                       21075, 21225, 37730,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 54650, 0, 3,
                                                                       52445, 36380, 52760,
                                                                       21525, 21735, 37955,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 55091, 0, 3,
                                                                       52760, 36605, 53075,
                                                                       21735, 21945, 38270,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 55532, 0, 3,
                                                                       53075, 36830, 53390,
                                                                       21945, 22155, 38585,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 55973, 0, 3,
                                                                       53390, 37055, 53705,
                                                                       22155, 22365, 38900,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 56414, 0, 3,
                                                                       53705, 37280, 54020,
                                                                       22365, 22575, 39215,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 56855, 0, 3,
                                                                       54020, 37505, 54335,
                                                                       22575, 22785, 39530,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 57296, 0, 3,
                                                                       54650, 37955, 55091,
                                                                       23205, 23485, 39845,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 57884, 0, 3,
                                                                       55091, 38270, 55532,
                                                                       23485, 23765, 40265,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 58472, 0, 3,
                                                                       55532, 38585, 55973,
                                                                       23765, 24045, 40685,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 59060, 0, 3,
                                                                       55973, 38900, 56414,
                                                                       24045, 24325, 41105,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 59648, 0, 3,
                                                                       56414, 39215, 56855,
                                                                       24325, 24605, 41525,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 60236, 0, 3,
                                                                       57296, 39845, 57884,
                                                                       25165, 25525, 41945,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 60992, 0, 3,
                                                                       57884, 40265, 58472,
                                                                       25525, 25885, 42485,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 61748, 0, 3,
                                                                       58472, 40685, 59060,
                                                                       25885, 26245, 43025,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 62504, 0, 3,
                                                                       59060, 41105, 59648,
                                                                       26245, 26605, 43565,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 63260, 0, 3,
                                                                       60236, 41945, 60992,
                                                                       27325, 27775, 44105,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 64205, 0, 3,
                                                                       60992, 42485, 61748,
                                                                       27775, 28225, 44780,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 65150, 0, 3,
                                                                       61748, 43025, 62504,
                                                                       28225, 28675, 45455,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 66095, 0, 3,
                                                                       63260, 44105, 64205,
                                                                       29575, 30125, 46130,
                                                                       ncols, gamma, p, q);

                    compute_prim_msh_three_center_electron_repulsion_0(buffer, 67250, 0, 3,
                                                                       64205, 44780, 65150,
                                                                       30125, 30675, 46955,
                                                                       ncols, gamma, p, q);

                    compute_prim_nsh_three_center_electron_repulsion_0(buffer, 68405, 0, 3,
                                                                       66095, 46130, 67250,
                                                                       31775, 32435, 47780,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 69791, 57296, 588, ncols);

                    simdfunc::contract_primitives(buffer, 70687, 60236, 756, ncols);

                    simdfunc::contract_primitives(buffer, 71839, 63260, 945, ncols);

                    simdfunc::contract_primitives(buffer, 73279, 66095, 1155, ncols);

                    simdfunc::contract_primitives(buffer, 75039, 68405, 1386, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 70379, 69791, 28, 1, nmax);

        simdtrf::transform_h_inner(buffer, 71443, 70687, 36, 1, nmax);

        simdtrf::transform_h_inner(buffer, 72784, 71839, 45, 1, nmax);

        simdtrf::transform_h_inner(buffer, 74434, 73279, 55, 1, nmax);

        simdtrf::transform_h_inner(buffer, 76425, 75039, 66, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 77151, 70379, 71443, 11,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 78075, 71443, 72784, 11,
                                             nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 79263, 72784, 74434, 11,
                                             nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 80748, 74434, 76425, 11,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 82563, 77151, 78075, 11,
                                             nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 84411, 78075, 79263, 11,
                                             nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 86787, 79263, 80748, 11,
                                             nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 89757, 82563, 84411, 11,
                                             nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 92837, 84411, 86787, 11,
                                             nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 96797, 89757, 92837, 11,
                                             nmax);

        simdtrf::transform_g_inner(buffer, 101417, 96797, 28, 11, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 101417, 99, nmax);
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
