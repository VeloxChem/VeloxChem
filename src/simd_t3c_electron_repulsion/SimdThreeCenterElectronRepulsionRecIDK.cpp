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


#include "SimdThreeCenterElectronRepulsionRecIDK.hpp"

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
#include "SimdTransferIP.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_idk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_idk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 98866, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 975 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 98866, 85807, 4884, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1492, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1495, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1498, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1501, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1504, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1507, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1510, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1513, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1516, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1519, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1522, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1525, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1528, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1531, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1534, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1537, 3, 7, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1546, 3, 8, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1555, 3, 9, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1564, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1573, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1582, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1591, 3, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1600, 3, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1609, 3, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1618, 3, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1627, 3, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1636, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1645, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1654, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1663, 3, 22, 64,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1681, 3, 25, 70,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1699, 3, 28, 76,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1717, 3, 31, 82,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1735, 3, 34, 88,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1753, 3, 37, 94,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1771, 3, 40, 100,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1789, 3, 43, 106,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1807, 3, 46, 112,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1825, 3, 49, 118,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1843, 3, 52, 124,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1861, 3, 55, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1879, 3, 58, 136,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1897, 3, 64, 142,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1927, 3, 70, 152,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1957, 3, 76, 162,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1987, 3, 82, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2017, 3, 88, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2047, 3, 94, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2077, 3, 100, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2107, 3, 106, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2137, 3, 112, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2167, 3, 118, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2197, 3, 124, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2227, 3, 130, 252,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2257, 3, 142, 262,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2302, 3, 152, 277,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2347, 3, 162, 292,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2392, 3, 172, 307,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2437, 3, 182, 322,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2482, 3, 192, 337,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2527, 3, 202, 352,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2572, 3, 212, 367,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2617, 3, 222, 382,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2662, 3, 232, 397,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2707, 3, 242, 412,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2752, 3, 262, 427,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2815, 3, 277, 448,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2878, 3, 292, 469,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2941, 3, 307, 490,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3004, 3, 322, 511,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3067, 3, 337, 532,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3130, 3, 352, 553,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3193, 3, 367, 574,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3256, 3, 382, 595,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3319, 3, 397, 616,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3382, 3, 427, 637,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3466, 3, 448, 665,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3550, 3, 469, 693,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3634, 3, 490, 721,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3718, 3, 511, 749,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3802, 3, 532, 777,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3886, 3, 553, 805,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3970, 3, 574, 833,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4054, 3, 595, 861,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4138, 3, 637, 889,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4246, 3, 665, 925,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4354, 3, 693, 961,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4462, 3, 721, 997,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4570, 3, 749,
                                                                       1033, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4678, 3, 777,
                                                                       1069, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4786, 3, 805,
                                                                       1105, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4894, 3, 833,
                                                                       1141, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5002, 3, 889,
                                                                       1177, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5137, 3, 925,
                                                                       1222, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5272, 3, 961,
                                                                       1267, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5407, 3, 997,
                                                                       1312, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5542, 3, 1033,
                                                                       1357, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5677, 3, 1069,
                                                                       1402, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 5812, 3, 1105,
                                                                       1447, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5947, 3, 7, 8,
                                                                       1498, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5953, 3, 8, 9,
                                                                       1501, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5959, 3, 9, 10,
                                                                       1504, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5965, 3, 10, 11,
                                                                       1507, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5971, 3, 11, 12,
                                                                       1510, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5977, 3, 12, 13,
                                                                       1513, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5983, 3, 13, 14,
                                                                       1516, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5989, 3, 14, 15,
                                                                       1519, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5995, 3, 15, 16,
                                                                       1522, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6001, 3, 16, 17,
                                                                       1525, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6007, 3, 17, 18,
                                                                       1528, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6013, 3, 18, 19,
                                                                       1531, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6019, 3, 19, 20,
                                                                       1534, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6025, 0, 3, 5947,
                                                                       1498, 5953, 1555, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6043, 0, 3, 5953,
                                                                       1501, 5959, 1564, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6061, 0, 3, 5959,
                                                                       1504, 5965, 1573, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6079, 0, 3, 5965,
                                                                       1507, 5971, 1582, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6097, 0, 3, 5971,
                                                                       1510, 5977, 1591, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6115, 0, 3, 5977,
                                                                       1513, 5983, 1600, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6133, 0, 3, 5983,
                                                                       1516, 5989, 1609, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6151, 0, 3, 5989,
                                                                       1519, 5995, 1618, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6169, 0, 3, 5995,
                                                                       1522, 6001, 1627, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6187, 0, 3, 6001,
                                                                       1525, 6007, 1636, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6205, 0, 3, 6007,
                                                                       1528, 6013, 1645, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6223, 0, 3, 6013,
                                                                       1531, 6019, 1654, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6241, 0, 3, 6025,
                                                                       1555, 6043, 64, 70, 1699,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6277, 0, 3, 6043,
                                                                       1564, 6061, 70, 76, 1717,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6313, 0, 3, 6061,
                                                                       1573, 6079, 76, 82, 1735,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6349, 0, 3, 6079,
                                                                       1582, 6097, 82, 88, 1753,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6385, 0, 3, 6097,
                                                                       1591, 6115, 88, 94, 1771,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6421, 0, 3, 6115,
                                                                       1600, 6133, 94, 100, 1789,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6457, 0, 3, 6133,
                                                                       1609, 6151, 100, 106,
                                                                       1807, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6493, 0, 3, 6151,
                                                                       1618, 6169, 106, 112,
                                                                       1825, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6529, 0, 3, 6169,
                                                                       1627, 6187, 112, 118,
                                                                       1843, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6565, 0, 3, 6187,
                                                                       1636, 6205, 118, 124,
                                                                       1861, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6601, 0, 3, 6205,
                                                                       1645, 6223, 124, 130,
                                                                       1879, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6637, 0, 3, 6241,
                                                                       1699, 6277, 142, 152,
                                                                       1957, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6697, 0, 3, 6277,
                                                                       1717, 6313, 152, 162,
                                                                       1987, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6757, 0, 3, 6313,
                                                                       1735, 6349, 162, 172,
                                                                       2017, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6817, 0, 3, 6349,
                                                                       1753, 6385, 172, 182,
                                                                       2047, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6877, 0, 3, 6385,
                                                                       1771, 6421, 182, 192,
                                                                       2077, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6937, 0, 3, 6421,
                                                                       1789, 6457, 192, 202,
                                                                       2107, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6997, 0, 3, 6457,
                                                                       1807, 6493, 202, 212,
                                                                       2137, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7057, 0, 3, 6493,
                                                                       1825, 6529, 212, 222,
                                                                       2167, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7117, 0, 3, 6529,
                                                                       1843, 6565, 222, 232,
                                                                       2197, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7177, 0, 3, 6565,
                                                                       1861, 6601, 232, 242,
                                                                       2227, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7237, 0, 3, 6637,
                                                                       1957, 6697, 262, 277,
                                                                       2347, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7327, 0, 3, 6697,
                                                                       1987, 6757, 277, 292,
                                                                       2392, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7417, 0, 3, 6757,
                                                                       2017, 6817, 292, 307,
                                                                       2437, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7507, 0, 3, 6817,
                                                                       2047, 6877, 307, 322,
                                                                       2482, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7597, 0, 3, 6877,
                                                                       2077, 6937, 322, 337,
                                                                       2527, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7687, 0, 3, 6937,
                                                                       2107, 6997, 337, 352,
                                                                       2572, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7777, 0, 3, 6997,
                                                                       2137, 7057, 352, 367,
                                                                       2617, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7867, 0, 3, 7057,
                                                                       2167, 7117, 367, 382,
                                                                       2662, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7957, 0, 3, 7117,
                                                                       2197, 7177, 382, 397,
                                                                       2707, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8047, 0, 3, 7237,
                                                                       2347, 7327, 427, 448,
                                                                       2878, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8173, 0, 3, 7327,
                                                                       2392, 7417, 448, 469,
                                                                       2941, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8299, 0, 3, 7417,
                                                                       2437, 7507, 469, 490,
                                                                       3004, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8425, 0, 3, 7507,
                                                                       2482, 7597, 490, 511,
                                                                       3067, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8551, 0, 3, 7597,
                                                                       2527, 7687, 511, 532,
                                                                       3130, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8677, 0, 3, 7687,
                                                                       2572, 7777, 532, 553,
                                                                       3193, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8803, 0, 3, 7777,
                                                                       2617, 7867, 553, 574,
                                                                       3256, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8929, 0, 3, 7867,
                                                                       2662, 7957, 574, 595,
                                                                       3319, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9055, 0, 3, 8047,
                                                                       2878, 8173, 637, 665,
                                                                       3550, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9223, 0, 3, 8173,
                                                                       2941, 8299, 665, 693,
                                                                       3634, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9391, 0, 3, 8299,
                                                                       3004, 8425, 693, 721,
                                                                       3718, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9559, 0, 3, 8425,
                                                                       3067, 8551, 721, 749,
                                                                       3802, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9727, 0, 3, 8551,
                                                                       3130, 8677, 749, 777,
                                                                       3886, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9895, 0, 3, 8677,
                                                                       3193, 8803, 777, 805,
                                                                       3970, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 10063, 0, 3, 8803,
                                                                       3256, 8929, 805, 833,
                                                                       4054, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 10231, 0, 3, 9055,
                                                                       3550, 9223, 889, 925,
                                                                       4354, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 10447, 0, 3, 9223,
                                                                       3634, 9391, 925, 961,
                                                                       4462, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 10663, 0, 3, 9391,
                                                                       3718, 9559, 961, 997,
                                                                       4570, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 10879, 0, 3, 9559,
                                                                       3802, 9727, 997, 1033,
                                                                       4678, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 11095, 0, 3, 9727,
                                                                       3886, 9895, 1033, 1069,
                                                                       4786, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 11311, 0, 3, 9895,
                                                                       3970, 10063, 1069, 1105,
                                                                       4894, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 11527, 0, 3,
                                                                       10231, 4354, 10447, 1177,
                                                                       1222, 5272, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 11797, 0, 3,
                                                                       10447, 4462, 10663, 1222,
                                                                       1267, 5407, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 12067, 0, 3,
                                                                       10663, 4570, 10879, 1267,
                                                                       1312, 5542, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 12337, 0, 3,
                                                                       10879, 4678, 11095, 1312,
                                                                       1357, 5677, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 12607, 0, 3,
                                                                       11095, 4786, 11311, 1357,
                                                                       1402, 5812, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12877, 3, 1492,
                                                                       1495, 5947, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12887, 3, 1495,
                                                                       1498, 5953, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12897, 3, 1498,
                                                                       1501, 5959, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12907, 3, 1501,
                                                                       1504, 5965, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12917, 3, 1504,
                                                                       1507, 5971, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12927, 3, 1507,
                                                                       1510, 5977, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12937, 3, 1510,
                                                                       1513, 5983, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12947, 3, 1513,
                                                                       1516, 5989, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12957, 3, 1516,
                                                                       1519, 5995, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12967, 3, 1519,
                                                                       1522, 6001, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12977, 3, 1522,
                                                                       1525, 6007, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12987, 3, 1525,
                                                                       1528, 6013, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12997, 3, 1528,
                                                                       1531, 6019, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13007, 0, 3,
                                                                       12877, 5947, 12887, 1537,
                                                                       1546, 6025, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13037, 0, 3,
                                                                       12887, 5953, 12897, 1546,
                                                                       1555, 6043, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13067, 0, 3,
                                                                       12897, 5959, 12907, 1555,
                                                                       1564, 6061, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13097, 0, 3,
                                                                       12907, 5965, 12917, 1564,
                                                                       1573, 6079, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13127, 0, 3,
                                                                       12917, 5971, 12927, 1573,
                                                                       1582, 6097, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13157, 0, 3,
                                                                       12927, 5977, 12937, 1582,
                                                                       1591, 6115, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13187, 0, 3,
                                                                       12937, 5983, 12947, 1591,
                                                                       1600, 6133, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13217, 0, 3,
                                                                       12947, 5989, 12957, 1600,
                                                                       1609, 6151, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13247, 0, 3,
                                                                       12957, 5995, 12967, 1609,
                                                                       1618, 6169, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13277, 0, 3,
                                                                       12967, 6001, 12977, 1618,
                                                                       1627, 6187, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13307, 0, 3,
                                                                       12977, 6007, 12987, 1627,
                                                                       1636, 6205, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13337, 0, 3,
                                                                       12987, 6013, 12997, 1636,
                                                                       1645, 6223, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13367, 0, 3,
                                                                       13007, 6025, 13037, 1663,
                                                                       1681, 6241, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13427, 0, 3,
                                                                       13037, 6043, 13067, 1681,
                                                                       1699, 6277, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13487, 0, 3,
                                                                       13067, 6061, 13097, 1699,
                                                                       1717, 6313, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13547, 0, 3,
                                                                       13097, 6079, 13127, 1717,
                                                                       1735, 6349, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13607, 0, 3,
                                                                       13127, 6097, 13157, 1735,
                                                                       1753, 6385, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13667, 0, 3,
                                                                       13157, 6115, 13187, 1753,
                                                                       1771, 6421, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13727, 0, 3,
                                                                       13187, 6133, 13217, 1771,
                                                                       1789, 6457, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13787, 0, 3,
                                                                       13217, 6151, 13247, 1789,
                                                                       1807, 6493, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13847, 0, 3,
                                                                       13247, 6169, 13277, 1807,
                                                                       1825, 6529, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13907, 0, 3,
                                                                       13277, 6187, 13307, 1825,
                                                                       1843, 6565, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13967, 0, 3,
                                                                       13307, 6205, 13337, 1843,
                                                                       1861, 6601, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14027, 0, 3,
                                                                       13367, 6241, 13427, 1897,
                                                                       1927, 6637, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14127, 0, 3,
                                                                       13427, 6277, 13487, 1927,
                                                                       1957, 6697, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14227, 0, 3,
                                                                       13487, 6313, 13547, 1957,
                                                                       1987, 6757, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14327, 0, 3,
                                                                       13547, 6349, 13607, 1987,
                                                                       2017, 6817, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14427, 0, 3,
                                                                       13607, 6385, 13667, 2017,
                                                                       2047, 6877, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14527, 0, 3,
                                                                       13667, 6421, 13727, 2047,
                                                                       2077, 6937, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14627, 0, 3,
                                                                       13727, 6457, 13787, 2077,
                                                                       2107, 6997, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14727, 0, 3,
                                                                       13787, 6493, 13847, 2107,
                                                                       2137, 7057, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14827, 0, 3,
                                                                       13847, 6529, 13907, 2137,
                                                                       2167, 7117, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14927, 0, 3,
                                                                       13907, 6565, 13967, 2167,
                                                                       2197, 7177, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15027, 0, 3,
                                                                       14027, 6637, 14127, 2257,
                                                                       2302, 7237, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15177, 0, 3,
                                                                       14127, 6697, 14227, 2302,
                                                                       2347, 7327, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15327, 0, 3,
                                                                       14227, 6757, 14327, 2347,
                                                                       2392, 7417, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15477, 0, 3,
                                                                       14327, 6817, 14427, 2392,
                                                                       2437, 7507, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15627, 0, 3,
                                                                       14427, 6877, 14527, 2437,
                                                                       2482, 7597, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15777, 0, 3,
                                                                       14527, 6937, 14627, 2482,
                                                                       2527, 7687, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 15927, 0, 3,
                                                                       14627, 6997, 14727, 2527,
                                                                       2572, 7777, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16077, 0, 3,
                                                                       14727, 7057, 14827, 2572,
                                                                       2617, 7867, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 16227, 0, 3,
                                                                       14827, 7117, 14927, 2617,
                                                                       2662, 7957, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 16377, 0, 3,
                                                                       15027, 7237, 15177, 2752,
                                                                       2815, 8047, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 16587, 0, 3,
                                                                       15177, 7327, 15327, 2815,
                                                                       2878, 8173, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 16797, 0, 3,
                                                                       15327, 7417, 15477, 2878,
                                                                       2941, 8299, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 17007, 0, 3,
                                                                       15477, 7507, 15627, 2941,
                                                                       3004, 8425, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 17217, 0, 3,
                                                                       15627, 7597, 15777, 3004,
                                                                       3067, 8551, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 17427, 0, 3,
                                                                       15777, 7687, 15927, 3067,
                                                                       3130, 8677, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 17637, 0, 3,
                                                                       15927, 7777, 16077, 3130,
                                                                       3193, 8803, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 17847, 0, 3,
                                                                       16077, 7867, 16227, 3193,
                                                                       3256, 8929, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 18057, 0, 3,
                                                                       16377, 8047, 16587, 3382,
                                                                       3466, 9055, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 18337, 0, 3,
                                                                       16587, 8173, 16797, 3466,
                                                                       3550, 9223, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 18617, 0, 3,
                                                                       16797, 8299, 17007, 3550,
                                                                       3634, 9391, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 18897, 0, 3,
                                                                       17007, 8425, 17217, 3634,
                                                                       3718, 9559, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 19177, 0, 3,
                                                                       17217, 8551, 17427, 3718,
                                                                       3802, 9727, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 19457, 0, 3,
                                                                       17427, 8677, 17637, 3802,
                                                                       3886, 9895, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 19737, 0, 3,
                                                                       17637, 8803, 17847, 3886,
                                                                       3970, 10063, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 20017, 0, 3,
                                                                       18057, 9055, 18337, 4138,
                                                                       4246, 10231, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 20377, 0, 3,
                                                                       18337, 9223, 18617, 4246,
                                                                       4354, 10447, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 20737, 0, 3,
                                                                       18617, 9391, 18897, 4354,
                                                                       4462, 10663, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 21097, 0, 3,
                                                                       18897, 9559, 19177, 4462,
                                                                       4570, 10879, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 21457, 0, 3,
                                                                       19177, 9727, 19457, 4570,
                                                                       4678, 11095, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 21817, 0, 3,
                                                                       19457, 9895, 19737, 4678,
                                                                       4786, 11311, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 22177, 0, 3,
                                                                       20017, 10231, 20377, 5002,
                                                                       5137, 11527, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 22627, 0, 3,
                                                                       20377, 10447, 20737, 5137,
                                                                       5272, 11797, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 23077, 0, 3,
                                                                       20737, 10663, 21097, 5272,
                                                                       5407, 12067, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 23527, 0, 3,
                                                                       21097, 10879, 21457, 5407,
                                                                       5542, 12337, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 23977, 0, 3,
                                                                       21457, 11095, 21817, 5542,
                                                                       5677, 12607, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24427, 3, 5947,
                                                                       5953, 12897, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24442, 3, 5953,
                                                                       5959, 12907, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24457, 3, 5959,
                                                                       5965, 12917, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24472, 3, 5965,
                                                                       5971, 12927, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24487, 3, 5971,
                                                                       5977, 12937, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24502, 3, 5977,
                                                                       5983, 12947, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24517, 3, 5983,
                                                                       5989, 12957, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24532, 3, 5989,
                                                                       5995, 12967, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24547, 3, 5995,
                                                                       6001, 12977, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24562, 3, 6001,
                                                                       6007, 12987, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 24577, 3, 6007,
                                                                       6013, 12997, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24592, 0, 3,
                                                                       24427, 12897, 24442, 6025,
                                                                       6043, 13067, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24637, 0, 3,
                                                                       24442, 12907, 24457, 6043,
                                                                       6061, 13097, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24682, 0, 3,
                                                                       24457, 12917, 24472, 6061,
                                                                       6079, 13127, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24727, 0, 3,
                                                                       24472, 12927, 24487, 6079,
                                                                       6097, 13157, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24772, 0, 3,
                                                                       24487, 12937, 24502, 6097,
                                                                       6115, 13187, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24817, 0, 3,
                                                                       24502, 12947, 24517, 6115,
                                                                       6133, 13217, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24862, 0, 3,
                                                                       24517, 12957, 24532, 6133,
                                                                       6151, 13247, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24907, 0, 3,
                                                                       24532, 12967, 24547, 6151,
                                                                       6169, 13277, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24952, 0, 3,
                                                                       24547, 12977, 24562, 6169,
                                                                       6187, 13307, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 24997, 0, 3,
                                                                       24562, 12987, 24577, 6187,
                                                                       6205, 13337, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25042, 0, 3,
                                                                       24592, 13067, 24637, 6241,
                                                                       6277, 13487, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25132, 0, 3,
                                                                       24637, 13097, 24682, 6277,
                                                                       6313, 13547, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25222, 0, 3,
                                                                       24682, 13127, 24727, 6313,
                                                                       6349, 13607, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25312, 0, 3,
                                                                       24727, 13157, 24772, 6349,
                                                                       6385, 13667, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25402, 0, 3,
                                                                       24772, 13187, 24817, 6385,
                                                                       6421, 13727, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25492, 0, 3,
                                                                       24817, 13217, 24862, 6421,
                                                                       6457, 13787, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25582, 0, 3,
                                                                       24862, 13247, 24907, 6457,
                                                                       6493, 13847, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25672, 0, 3,
                                                                       24907, 13277, 24952, 6493,
                                                                       6529, 13907, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 25762, 0, 3,
                                                                       24952, 13307, 24997, 6529,
                                                                       6565, 13967, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 25852, 0, 3,
                                                                       25042, 13487, 25132, 6637,
                                                                       6697, 14227, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26002, 0, 3,
                                                                       25132, 13547, 25222, 6697,
                                                                       6757, 14327, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26152, 0, 3,
                                                                       25222, 13607, 25312, 6757,
                                                                       6817, 14427, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26302, 0, 3,
                                                                       25312, 13667, 25402, 6817,
                                                                       6877, 14527, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26452, 0, 3,
                                                                       25402, 13727, 25492, 6877,
                                                                       6937, 14627, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26602, 0, 3,
                                                                       25492, 13787, 25582, 6937,
                                                                       6997, 14727, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26752, 0, 3,
                                                                       25582, 13847, 25672, 6997,
                                                                       7057, 14827, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 26902, 0, 3,
                                                                       25672, 13907, 25762, 7057,
                                                                       7117, 14927, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27052, 0, 3,
                                                                       25852, 14227, 26002, 7237,
                                                                       7327, 15327, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27277, 0, 3,
                                                                       26002, 14327, 26152, 7327,
                                                                       7417, 15477, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27502, 0, 3,
                                                                       26152, 14427, 26302, 7417,
                                                                       7507, 15627, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27727, 0, 3,
                                                                       26302, 14527, 26452, 7507,
                                                                       7597, 15777, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 27952, 0, 3,
                                                                       26452, 14627, 26602, 7597,
                                                                       7687, 15927, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 28177, 0, 3,
                                                                       26602, 14727, 26752, 7687,
                                                                       7777, 16077, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 28402, 0, 3,
                                                                       26752, 14827, 26902, 7777,
                                                                       7867, 16227, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 28627, 0, 3,
                                                                       27052, 15327, 27277, 8047,
                                                                       8173, 16797, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 28942, 0, 3,
                                                                       27277, 15477, 27502, 8173,
                                                                       8299, 17007, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 29257, 0, 3,
                                                                       27502, 15627, 27727, 8299,
                                                                       8425, 17217, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 29572, 0, 3,
                                                                       27727, 15777, 27952, 8425,
                                                                       8551, 17427, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 29887, 0, 3,
                                                                       27952, 15927, 28177, 8551,
                                                                       8677, 17637, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 30202, 0, 3,
                                                                       28177, 16077, 28402, 8677,
                                                                       8803, 17847, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 30517, 0, 3,
                                                                       28627, 16797, 28942, 9055,
                                                                       9223, 18617, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 30937, 0, 3,
                                                                       28942, 17007, 29257, 9223,
                                                                       9391, 18897, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 31357, 0, 3,
                                                                       29257, 17217, 29572, 9391,
                                                                       9559, 19177, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 31777, 0, 3,
                                                                       29572, 17427, 29887, 9559,
                                                                       9727, 19457, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 32197, 0, 3,
                                                                       29887, 17637, 30202, 9727,
                                                                       9895, 19737, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 32617, 0, 3,
                                                                       30517, 18617, 30937,
                                                                       10231, 10447, 20737,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 33157, 0, 3,
                                                                       30937, 18897, 31357,
                                                                       10447, 10663, 21097,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 33697, 0, 3,
                                                                       31357, 19177, 31777,
                                                                       10663, 10879, 21457,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 34237, 0, 3,
                                                                       31777, 19457, 32197,
                                                                       10879, 11095, 21817,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 34777, 0, 3,
                                                                       32617, 20737, 33157,
                                                                       11527, 11797, 23077,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 35452, 0, 3,
                                                                       33157, 21097, 33697,
                                                                       11797, 12067, 23527,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 36127, 0, 3,
                                                                       33697, 21457, 34237,
                                                                       12067, 12337, 23977,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36802, 3, 12877,
                                                                       12887, 24427, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36823, 3, 12887,
                                                                       12897, 24442, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36844, 3, 12897,
                                                                       12907, 24457, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36865, 3, 12907,
                                                                       12917, 24472, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36886, 3, 12917,
                                                                       12927, 24487, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36907, 3, 12927,
                                                                       12937, 24502, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36928, 3, 12937,
                                                                       12947, 24517, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36949, 3, 12947,
                                                                       12957, 24532, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36970, 3, 12957,
                                                                       12967, 24547, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 36991, 3, 12967,
                                                                       12977, 24562, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 37012, 3, 12977,
                                                                       12987, 24577, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37033, 0, 3,
                                                                       36802, 24427, 36823,
                                                                       13007, 13037, 24592,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37096, 0, 3,
                                                                       36823, 24442, 36844,
                                                                       13037, 13067, 24637,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37159, 0, 3,
                                                                       36844, 24457, 36865,
                                                                       13067, 13097, 24682,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37222, 0, 3,
                                                                       36865, 24472, 36886,
                                                                       13097, 13127, 24727,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37285, 0, 3,
                                                                       36886, 24487, 36907,
                                                                       13127, 13157, 24772,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37348, 0, 3,
                                                                       36907, 24502, 36928,
                                                                       13157, 13187, 24817,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37411, 0, 3,
                                                                       36928, 24517, 36949,
                                                                       13187, 13217, 24862,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37474, 0, 3,
                                                                       36949, 24532, 36970,
                                                                       13217, 13247, 24907,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37537, 0, 3,
                                                                       36970, 24547, 36991,
                                                                       13247, 13277, 24952,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 37600, 0, 3,
                                                                       36991, 24562, 37012,
                                                                       13277, 13307, 24997,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 37663, 0, 3,
                                                                       37033, 24592, 37096,
                                                                       13367, 13427, 25042,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 37789, 0, 3,
                                                                       37096, 24637, 37159,
                                                                       13427, 13487, 25132,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 37915, 0, 3,
                                                                       37159, 24682, 37222,
                                                                       13487, 13547, 25222,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38041, 0, 3,
                                                                       37222, 24727, 37285,
                                                                       13547, 13607, 25312,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38167, 0, 3,
                                                                       37285, 24772, 37348,
                                                                       13607, 13667, 25402,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38293, 0, 3,
                                                                       37348, 24817, 37411,
                                                                       13667, 13727, 25492,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38419, 0, 3,
                                                                       37411, 24862, 37474,
                                                                       13727, 13787, 25582,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38545, 0, 3,
                                                                       37474, 24907, 37537,
                                                                       13787, 13847, 25672,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38671, 0, 3,
                                                                       37537, 24952, 37600,
                                                                       13847, 13907, 25762,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 38797, 0, 3,
                                                                       37663, 25042, 37789,
                                                                       14027, 14127, 25852,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39007, 0, 3,
                                                                       37789, 25132, 37915,
                                                                       14127, 14227, 26002,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39217, 0, 3,
                                                                       37915, 25222, 38041,
                                                                       14227, 14327, 26152,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39427, 0, 3,
                                                                       38041, 25312, 38167,
                                                                       14327, 14427, 26302,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39637, 0, 3,
                                                                       38167, 25402, 38293,
                                                                       14427, 14527, 26452,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 39847, 0, 3,
                                                                       38293, 25492, 38419,
                                                                       14527, 14627, 26602,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 40057, 0, 3,
                                                                       38419, 25582, 38545,
                                                                       14627, 14727, 26752,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 40267, 0, 3,
                                                                       38545, 25672, 38671,
                                                                       14727, 14827, 26902,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 40477, 0, 3,
                                                                       38797, 25852, 39007,
                                                                       15027, 15177, 27052,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 40792, 0, 3,
                                                                       39007, 26002, 39217,
                                                                       15177, 15327, 27277,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 41107, 0, 3,
                                                                       39217, 26152, 39427,
                                                                       15327, 15477, 27502,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 41422, 0, 3,
                                                                       39427, 26302, 39637,
                                                                       15477, 15627, 27727,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 41737, 0, 3,
                                                                       39637, 26452, 39847,
                                                                       15627, 15777, 27952,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 42052, 0, 3,
                                                                       39847, 26602, 40057,
                                                                       15777, 15927, 28177,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 42367, 0, 3,
                                                                       40057, 26752, 40267,
                                                                       15927, 16077, 28402,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 42682, 0, 3,
                                                                       40477, 27052, 40792,
                                                                       16377, 16587, 28627,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 43123, 0, 3,
                                                                       40792, 27277, 41107,
                                                                       16587, 16797, 28942,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 43564, 0, 3,
                                                                       41107, 27502, 41422,
                                                                       16797, 17007, 29257,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 44005, 0, 3,
                                                                       41422, 27727, 41737,
                                                                       17007, 17217, 29572,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 44446, 0, 3,
                                                                       41737, 27952, 42052,
                                                                       17217, 17427, 29887,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 44887, 0, 3,
                                                                       42052, 28177, 42367,
                                                                       17427, 17637, 30202,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 45328, 0, 3,
                                                                       42682, 28627, 43123,
                                                                       18057, 18337, 30517,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 45916, 0, 3,
                                                                       43123, 28942, 43564,
                                                                       18337, 18617, 30937,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 46504, 0, 3,
                                                                       43564, 29257, 44005,
                                                                       18617, 18897, 31357,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 47092, 0, 3,
                                                                       44005, 29572, 44446,
                                                                       18897, 19177, 31777,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 47680, 0, 3,
                                                                       44446, 29887, 44887,
                                                                       19177, 19457, 32197,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 48268, 0, 3,
                                                                       45328, 30517, 45916,
                                                                       20017, 20377, 32617,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 49024, 0, 3,
                                                                       45916, 30937, 46504,
                                                                       20377, 20737, 33157,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 49780, 0, 3,
                                                                       46504, 31357, 47092,
                                                                       20737, 21097, 33697,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 50536, 0, 3,
                                                                       47092, 31777, 47680,
                                                                       21097, 21457, 34237,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 51292, 0, 3,
                                                                       48268, 32617, 49024,
                                                                       22177, 22627, 34777,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 52237, 0, 3,
                                                                       49024, 33157, 49780,
                                                                       22627, 23077, 35452,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 53182, 0, 3,
                                                                       49780, 33697, 50536,
                                                                       23077, 23527, 36127,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54127, 3, 24427,
                                                                       24442, 36844, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54155, 3, 24442,
                                                                       24457, 36865, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54183, 3, 24457,
                                                                       24472, 36886, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54211, 3, 24472,
                                                                       24487, 36907, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54239, 3, 24487,
                                                                       24502, 36928, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54267, 3, 24502,
                                                                       24517, 36949, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54295, 3, 24517,
                                                                       24532, 36970, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54323, 3, 24532,
                                                                       24547, 36991, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 54351, 3, 24547,
                                                                       24562, 37012, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54379, 0, 3,
                                                                       54127, 36844, 54155,
                                                                       24592, 24637, 37159,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54463, 0, 3,
                                                                       54155, 36865, 54183,
                                                                       24637, 24682, 37222,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54547, 0, 3,
                                                                       54183, 36886, 54211,
                                                                       24682, 24727, 37285,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54631, 0, 3,
                                                                       54211, 36907, 54239,
                                                                       24727, 24772, 37348,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54715, 0, 3,
                                                                       54239, 36928, 54267,
                                                                       24772, 24817, 37411,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54799, 0, 3,
                                                                       54267, 36949, 54295,
                                                                       24817, 24862, 37474,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54883, 0, 3,
                                                                       54295, 36970, 54323,
                                                                       24862, 24907, 37537,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 54967, 0, 3,
                                                                       54323, 36991, 54351,
                                                                       24907, 24952, 37600,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 55051, 0, 3,
                                                                       54379, 37159, 54463,
                                                                       25042, 25132, 37915,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 55219, 0, 3,
                                                                       54463, 37222, 54547,
                                                                       25132, 25222, 38041,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 55387, 0, 3,
                                                                       54547, 37285, 54631,
                                                                       25222, 25312, 38167,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 55555, 0, 3,
                                                                       54631, 37348, 54715,
                                                                       25312, 25402, 38293,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 55723, 0, 3,
                                                                       54715, 37411, 54799,
                                                                       25402, 25492, 38419,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 55891, 0, 3,
                                                                       54799, 37474, 54883,
                                                                       25492, 25582, 38545,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 56059, 0, 3,
                                                                       54883, 37537, 54967,
                                                                       25582, 25672, 38671,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 56227, 0, 3,
                                                                       55051, 37915, 55219,
                                                                       25852, 26002, 39217,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 56507, 0, 3,
                                                                       55219, 38041, 55387,
                                                                       26002, 26152, 39427,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 56787, 0, 3,
                                                                       55387, 38167, 55555,
                                                                       26152, 26302, 39637,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 57067, 0, 3,
                                                                       55555, 38293, 55723,
                                                                       26302, 26452, 39847,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 57347, 0, 3,
                                                                       55723, 38419, 55891,
                                                                       26452, 26602, 40057,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 57627, 0, 3,
                                                                       55891, 38545, 56059,
                                                                       26602, 26752, 40267,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 57907, 0, 3,
                                                                       56227, 39217, 56507,
                                                                       27052, 27277, 41107,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 58327, 0, 3,
                                                                       56507, 39427, 56787,
                                                                       27277, 27502, 41422,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 58747, 0, 3,
                                                                       56787, 39637, 57067,
                                                                       27502, 27727, 41737,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 59167, 0, 3,
                                                                       57067, 39847, 57347,
                                                                       27727, 27952, 42052,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 59587, 0, 3,
                                                                       57347, 40057, 57627,
                                                                       27952, 28177, 42367,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 60007, 0, 3,
                                                                       57907, 41107, 58327,
                                                                       28627, 28942, 43564,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 60595, 0, 3,
                                                                       58327, 41422, 58747,
                                                                       28942, 29257, 44005,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 61183, 0, 3,
                                                                       58747, 41737, 59167,
                                                                       29257, 29572, 44446,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 61771, 0, 3,
                                                                       59167, 42052, 59587,
                                                                       29572, 29887, 44887,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 62359, 0, 3,
                                                                       60007, 43564, 60595,
                                                                       30517, 30937, 46504,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 63143, 0, 3,
                                                                       60595, 44005, 61183,
                                                                       30937, 31357, 47092,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 63927, 0, 3,
                                                                       61183, 44446, 61771,
                                                                       31357, 31777, 47680,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 64711, 0, 3,
                                                                       62359, 46504, 63143,
                                                                       32617, 33157, 49780,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 65719, 0, 3,
                                                                       63143, 47092, 63927,
                                                                       33157, 33697, 50536,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 66727, 0, 3,
                                                                       64711, 49780, 65719,
                                                                       34777, 35452, 53182,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 67987, 3, 36802,
                                                                       36823, 54127, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68023, 3, 36823,
                                                                       36844, 54155, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68059, 3, 36844,
                                                                       36865, 54183, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68095, 3, 36865,
                                                                       36886, 54211, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68131, 3, 36886,
                                                                       36907, 54239, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68167, 3, 36907,
                                                                       36928, 54267, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68203, 3, 36928,
                                                                       36949, 54295, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68239, 3, 36949,
                                                                       36970, 54323, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 68275, 3, 36970,
                                                                       36991, 54351, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68311, 0, 3,
                                                                       67987, 54127, 68023,
                                                                       37033, 37096, 54379,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68419, 0, 3,
                                                                       68023, 54155, 68059,
                                                                       37096, 37159, 54463,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68527, 0, 3,
                                                                       68059, 54183, 68095,
                                                                       37159, 37222, 54547,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68635, 0, 3,
                                                                       68095, 54211, 68131,
                                                                       37222, 37285, 54631,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68743, 0, 3,
                                                                       68131, 54239, 68167,
                                                                       37285, 37348, 54715,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68851, 0, 3,
                                                                       68167, 54267, 68203,
                                                                       37348, 37411, 54799,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 68959, 0, 3,
                                                                       68203, 54295, 68239,
                                                                       37411, 37474, 54883,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 69067, 0, 3,
                                                                       68239, 54323, 68275,
                                                                       37474, 37537, 54967,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 69175, 0, 3,
                                                                       68311, 54379, 68419,
                                                                       37663, 37789, 55051,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 69391, 0, 3,
                                                                       68419, 54463, 68527,
                                                                       37789, 37915, 55219,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 69607, 0, 3,
                                                                       68527, 54547, 68635,
                                                                       37915, 38041, 55387,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 69823, 0, 3,
                                                                       68635, 54631, 68743,
                                                                       38041, 38167, 55555,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 70039, 0, 3,
                                                                       68743, 54715, 68851,
                                                                       38167, 38293, 55723,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 70255, 0, 3,
                                                                       68851, 54799, 68959,
                                                                       38293, 38419, 55891,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 70471, 0, 3,
                                                                       68959, 54883, 69067,
                                                                       38419, 38545, 56059,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 70687, 0, 3,
                                                                       69175, 55051, 69391,
                                                                       38797, 39007, 56227,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 71047, 0, 3,
                                                                       69391, 55219, 69607,
                                                                       39007, 39217, 56507,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 71407, 0, 3,
                                                                       69607, 55387, 69823,
                                                                       39217, 39427, 56787,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 71767, 0, 3,
                                                                       69823, 55555, 70039,
                                                                       39427, 39637, 57067,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 72127, 0, 3,
                                                                       70039, 55723, 70255,
                                                                       39637, 39847, 57347,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 72487, 0, 3,
                                                                       70255, 55891, 70471,
                                                                       39847, 40057, 57627,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 72847, 0, 3,
                                                                       70687, 56227, 71047,
                                                                       40477, 40792, 57907,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 73387, 0, 3,
                                                                       71047, 56507, 71407,
                                                                       40792, 41107, 58327,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 73927, 0, 3,
                                                                       71407, 56787, 71767,
                                                                       41107, 41422, 58747,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 74467, 0, 3,
                                                                       71767, 57067, 72127,
                                                                       41422, 41737, 59167,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 75007, 0, 3,
                                                                       72127, 57347, 72487,
                                                                       41737, 42052, 59587,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 75547, 0, 3,
                                                                       72847, 57907, 73387,
                                                                       42682, 43123, 60007,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 76303, 0, 3,
                                                                       73387, 58327, 73927,
                                                                       43123, 43564, 60595,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 77059, 0, 3,
                                                                       73927, 58747, 74467,
                                                                       43564, 44005, 61183,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 77815, 0, 3,
                                                                       74467, 59167, 75007,
                                                                       44005, 44446, 61771,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 78571, 0, 3,
                                                                       75547, 60007, 76303,
                                                                       45328, 45916, 62359,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 79579, 0, 3,
                                                                       76303, 60595, 77059,
                                                                       45916, 46504, 63143,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 80587, 0, 3,
                                                                       77059, 61183, 77815,
                                                                       46504, 47092, 63927,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 81595, 0, 3,
                                                                       78571, 62359, 79579,
                                                                       48268, 49024, 64711,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksk_three_center_electron_repulsion_0(buffer, 82891, 0, 3,
                                                                       79579, 63143, 80587,
                                                                       49024, 49780, 65719,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsk_three_center_electron_repulsion_0(buffer, 84187, 0, 3,
                                                                       81595, 64711, 82891,
                                                                       51292, 52237, 66727,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 85807, 78571, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 87235, 81595, 1296, ncols);

                    simdfunc::contract_primitives(buffer, 89071, 84187, 1620, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 86815, 85807, 28, 1, nmax);

        simdtrf::transform_k_inner(buffer, 88531, 87235, 36, 1, nmax);

        simdtrf::transform_k_inner(buffer, 90691, 89071, 45, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 91366, 86815, 88531, 15,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 92626, 88531, 90691, 15,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 94246, 91366, 92626, 15,
                                             nmax);

        simdtrf::transform_d_inner(buffer, 96766, 94246, 28, 15, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 96766, 75, nmax);
    }

    for (size_t m = 0; m < 975; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
