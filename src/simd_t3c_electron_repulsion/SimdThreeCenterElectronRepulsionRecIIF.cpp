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


#include "SimdThreeCenterElectronRepulsionRecIIF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecNSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecOSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecQSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIF.hpp"
#include "SimdTransferIG.hpp"
#include "SimdTransferIH.hpp"
#include "SimdTransferII.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKD.hpp"
#include "SimdTransferKF.hpp"
#include "SimdTransferKG.hpp"
#include "SimdTransferKH.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransferLD.hpp"
#include "SimdTransferLF.hpp"
#include "SimdTransferLG.hpp"
#include "SimdTransferLP.hpp"
#include "SimdTransferMD.hpp"
#include "SimdTransferMF.hpp"
#include "SimdTransferMP.hpp"
#include "SimdTransferND.hpp"
#include "SimdTransferNP.hpp"
#include "SimdTransferOP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_iif_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_iif_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 103327, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1183 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 103327, 40047, 6146, dimensions);

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

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2152, 0, 3, 1492,
                                                                       1547, 1822, 1888, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2230, 0, 3, 1547,
                                                                       1602, 1888, 1954, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2308, 0, 3, 1602,
                                                                       1657, 1954, 2020, ncols,
                                                                       gamma, p, q);

                    compute_prim_oss_three_center_electron_repulsion_0(buffer, 2386, 0, 3, 1657,
                                                                       1712, 2020, 2086, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 2464, 0, 3, 1822,
                                                                       1888, 2152, 2230, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 2555, 0, 3, 1888,
                                                                       1954, 2230, 2308, ncols,
                                                                       gamma, p, q);

                    compute_prim_qss_three_center_electron_repulsion_0(buffer, 2646, 0, 3, 1954,
                                                                       2020, 2308, 2386, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2737, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2740, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2743, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2746, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2749, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2752, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2755, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2758, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2761, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2764, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2767, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2770, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2773, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2776, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2779, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2782, 3, 7, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2791, 3, 8, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2800, 3, 9, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2809, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2818, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2827, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2836, 3, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2845, 3, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2854, 3, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2863, 3, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2872, 3, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2881, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2890, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2899, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2908, 3, 22, 64,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2926, 3, 25, 70,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2944, 3, 28, 76,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2962, 3, 31, 82,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2980, 3, 34, 88,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2998, 3, 37, 94,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3016, 3, 40, 100,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3034, 3, 43, 106,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3052, 3, 46, 112,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3070, 3, 49, 118,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3088, 3, 52, 124,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3106, 3, 55, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3124, 3, 58, 136,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3142, 3, 64, 142,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3172, 3, 70, 152,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3202, 3, 76, 162,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3232, 3, 82, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3262, 3, 88, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3292, 3, 94, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3322, 3, 100, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3352, 3, 106, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3382, 3, 112, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3412, 3, 118, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3442, 3, 124, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3472, 3, 130, 252,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3502, 3, 142, 262,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3547, 3, 152, 277,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3592, 3, 162, 292,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3637, 3, 172, 307,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3682, 3, 182, 322,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3727, 3, 192, 337,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3772, 3, 202, 352,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3817, 3, 212, 367,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3862, 3, 222, 382,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3907, 3, 232, 397,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3952, 3, 242, 412,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3997, 3, 262, 427,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4060, 3, 277, 448,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4123, 3, 292, 469,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4186, 3, 307, 490,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4249, 3, 322, 511,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4312, 3, 337, 532,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4375, 3, 352, 553,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4438, 3, 367, 574,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4501, 3, 382, 595,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4564, 3, 397, 616,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4627, 3, 427, 637,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4711, 3, 448, 665,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4795, 3, 469, 693,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4879, 3, 490, 721,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4963, 3, 511, 749,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5047, 3, 532, 777,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5131, 3, 553, 805,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5215, 3, 574, 833,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5299, 3, 595, 861,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5383, 3, 637, 889,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5491, 3, 665, 925,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5599, 3, 693, 961,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5707, 3, 721, 997,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5815, 3, 749,
                                                                       1033, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5923, 3, 777,
                                                                       1069, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6031, 3, 805,
                                                                       1105, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6139, 3, 833,
                                                                       1141, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6247, 3, 889,
                                                                       1177, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6382, 3, 925,
                                                                       1222, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6517, 3, 961,
                                                                       1267, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6652, 3, 997,
                                                                       1312, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6787, 3, 1033,
                                                                       1357, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 6922, 3, 1069,
                                                                       1402, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 7057, 3, 1105,
                                                                       1447, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7192, 3, 1177,
                                                                       1492, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7357, 3, 1222,
                                                                       1547, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7522, 3, 1267,
                                                                       1602, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7687, 3, 1312,
                                                                       1657, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 7852, 3, 1357,
                                                                       1712, ncols, p, q);

                    compute_prim_msp_three_center_electron_repulsion_0(buffer, 8017, 3, 1402,
                                                                       1767, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8182, 3, 1492,
                                                                       1822, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8380, 3, 1547,
                                                                       1888, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8578, 3, 1602,
                                                                       1954, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8776, 3, 1657,
                                                                       2020, ncols, p, q);

                    compute_prim_nsp_three_center_electron_repulsion_0(buffer, 8974, 3, 1712,
                                                                       2086, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 9172, 3, 1822,
                                                                       2152, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 9406, 3, 1888,
                                                                       2230, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 9640, 3, 1954,
                                                                       2308, ncols, p, q);

                    compute_prim_osp_three_center_electron_repulsion_0(buffer, 9874, 3, 2020,
                                                                       2386, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 10108, 3, 2152,
                                                                       2464, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 10381, 3, 2230,
                                                                       2555, ncols, p, q);

                    compute_prim_qsp_three_center_electron_repulsion_0(buffer, 10654, 3, 2308,
                                                                       2646, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10927, 3, 7, 8,
                                                                       2743, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10933, 3, 8, 9,
                                                                       2746, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10939, 3, 9, 10,
                                                                       2749, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10945, 3, 10, 11,
                                                                       2752, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10951, 3, 11, 12,
                                                                       2755, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10957, 3, 12, 13,
                                                                       2758, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10963, 3, 13, 14,
                                                                       2761, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10969, 3, 14, 15,
                                                                       2764, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10975, 3, 15, 16,
                                                                       2767, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10981, 3, 16, 17,
                                                                       2770, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10987, 3, 17, 18,
                                                                       2773, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10993, 3, 18, 19,
                                                                       2776, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 10999, 3, 19, 20,
                                                                       2779, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11005, 0, 3,
                                                                       10927, 2743, 10933, 2800,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11023, 0, 3,
                                                                       10933, 2746, 10939, 2809,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11041, 0, 3,
                                                                       10939, 2749, 10945, 2818,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11059, 0, 3,
                                                                       10945, 2752, 10951, 2827,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11077, 0, 3,
                                                                       10951, 2755, 10957, 2836,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11095, 0, 3,
                                                                       10957, 2758, 10963, 2845,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11113, 0, 3,
                                                                       10963, 2761, 10969, 2854,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11131, 0, 3,
                                                                       10969, 2764, 10975, 2863,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11149, 0, 3,
                                                                       10975, 2767, 10981, 2872,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11167, 0, 3,
                                                                       10981, 2770, 10987, 2881,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11185, 0, 3,
                                                                       10987, 2773, 10993, 2890,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 11203, 0, 3,
                                                                       10993, 2776, 10999, 2899,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11221, 0, 3,
                                                                       11005, 2800, 11023, 64,
                                                                       70, 2944, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11257, 0, 3,
                                                                       11023, 2809, 11041, 70,
                                                                       76, 2962, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11293, 0, 3,
                                                                       11041, 2818, 11059, 76,
                                                                       82, 2980, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11329, 0, 3,
                                                                       11059, 2827, 11077, 82,
                                                                       88, 2998, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11365, 0, 3,
                                                                       11077, 2836, 11095, 88,
                                                                       94, 3016, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11401, 0, 3,
                                                                       11095, 2845, 11113, 94,
                                                                       100, 3034, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11437, 0, 3,
                                                                       11113, 2854, 11131, 100,
                                                                       106, 3052, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11473, 0, 3,
                                                                       11131, 2863, 11149, 106,
                                                                       112, 3070, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11509, 0, 3,
                                                                       11149, 2872, 11167, 112,
                                                                       118, 3088, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11545, 0, 3,
                                                                       11167, 2881, 11185, 118,
                                                                       124, 3106, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 11581, 0, 3,
                                                                       11185, 2890, 11203, 124,
                                                                       130, 3124, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11617, 0, 3,
                                                                       11221, 2944, 11257, 142,
                                                                       152, 3202, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11677, 0, 3,
                                                                       11257, 2962, 11293, 152,
                                                                       162, 3232, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11737, 0, 3,
                                                                       11293, 2980, 11329, 162,
                                                                       172, 3262, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11797, 0, 3,
                                                                       11329, 2998, 11365, 172,
                                                                       182, 3292, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11857, 0, 3,
                                                                       11365, 3016, 11401, 182,
                                                                       192, 3322, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11917, 0, 3,
                                                                       11401, 3034, 11437, 192,
                                                                       202, 3352, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11977, 0, 3,
                                                                       11437, 3052, 11473, 202,
                                                                       212, 3382, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12037, 0, 3,
                                                                       11473, 3070, 11509, 212,
                                                                       222, 3412, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12097, 0, 3,
                                                                       11509, 3088, 11545, 222,
                                                                       232, 3442, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12157, 0, 3,
                                                                       11545, 3106, 11581, 232,
                                                                       242, 3472, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12217, 0, 3,
                                                                       11617, 3202, 11677, 262,
                                                                       277, 3592, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12307, 0, 3,
                                                                       11677, 3232, 11737, 277,
                                                                       292, 3637, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12397, 0, 3,
                                                                       11737, 3262, 11797, 292,
                                                                       307, 3682, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12487, 0, 3,
                                                                       11797, 3292, 11857, 307,
                                                                       322, 3727, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12577, 0, 3,
                                                                       11857, 3322, 11917, 322,
                                                                       337, 3772, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12667, 0, 3,
                                                                       11917, 3352, 11977, 337,
                                                                       352, 3817, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12757, 0, 3,
                                                                       11977, 3382, 12037, 352,
                                                                       367, 3862, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12847, 0, 3,
                                                                       12037, 3412, 12097, 367,
                                                                       382, 3907, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12937, 0, 3,
                                                                       12097, 3442, 12157, 382,
                                                                       397, 3952, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13027, 0, 3,
                                                                       12217, 3592, 12307, 427,
                                                                       448, 4123, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13153, 0, 3,
                                                                       12307, 3637, 12397, 448,
                                                                       469, 4186, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13279, 0, 3,
                                                                       12397, 3682, 12487, 469,
                                                                       490, 4249, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13405, 0, 3,
                                                                       12487, 3727, 12577, 490,
                                                                       511, 4312, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13531, 0, 3,
                                                                       12577, 3772, 12667, 511,
                                                                       532, 4375, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13657, 0, 3,
                                                                       12667, 3817, 12757, 532,
                                                                       553, 4438, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13783, 0, 3,
                                                                       12757, 3862, 12847, 553,
                                                                       574, 4501, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 13909, 0, 3,
                                                                       12847, 3907, 12937, 574,
                                                                       595, 4564, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14035, 0, 3,
                                                                       13027, 4123, 13153, 637,
                                                                       665, 4795, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14203, 0, 3,
                                                                       13153, 4186, 13279, 665,
                                                                       693, 4879, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14371, 0, 3,
                                                                       13279, 4249, 13405, 693,
                                                                       721, 4963, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14539, 0, 3,
                                                                       13405, 4312, 13531, 721,
                                                                       749, 5047, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14707, 0, 3,
                                                                       13531, 4375, 13657, 749,
                                                                       777, 5131, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14875, 0, 3,
                                                                       13657, 4438, 13783, 777,
                                                                       805, 5215, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 15043, 0, 3,
                                                                       13783, 4501, 13909, 805,
                                                                       833, 5299, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15211, 0, 3,
                                                                       14035, 4795, 14203, 889,
                                                                       925, 5599, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15427, 0, 3,
                                                                       14203, 4879, 14371, 925,
                                                                       961, 5707, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15643, 0, 3,
                                                                       14371, 4963, 14539, 961,
                                                                       997, 5815, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15859, 0, 3,
                                                                       14539, 5047, 14707, 997,
                                                                       1033, 5923, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16075, 0, 3,
                                                                       14707, 5131, 14875, 1033,
                                                                       1069, 6031, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16291, 0, 3,
                                                                       14875, 5215, 15043, 1069,
                                                                       1105, 6139, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16507, 0, 3,
                                                                       15211, 5599, 15427, 1177,
                                                                       1222, 6517, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 16777, 0, 3,
                                                                       15427, 5707, 15643, 1222,
                                                                       1267, 6652, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17047, 0, 3,
                                                                       15643, 5815, 15859, 1267,
                                                                       1312, 6787, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17317, 0, 3,
                                                                       15859, 5923, 16075, 1312,
                                                                       1357, 6922, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 17587, 0, 3,
                                                                       16075, 6031, 16291, 1357,
                                                                       1402, 7057, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 17857, 0, 3,
                                                                       16507, 6517, 16777, 1492,
                                                                       1547, 7522, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 18187, 0, 3,
                                                                       16777, 6652, 17047, 1547,
                                                                       1602, 7687, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 18517, 0, 3,
                                                                       17047, 6787, 17317, 1602,
                                                                       1657, 7852, ncols, gamma,
                                                                       p, q);

                    compute_prim_msd_three_center_electron_repulsion_0(buffer, 18847, 0, 3,
                                                                       17317, 6922, 17587, 1657,
                                                                       1712, 8017, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 19177, 0, 3,
                                                                       17857, 7522, 18187, 1822,
                                                                       1888, 8578, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 19573, 0, 3,
                                                                       18187, 7687, 18517, 1888,
                                                                       1954, 8776, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsd_three_center_electron_repulsion_0(buffer, 19969, 0, 3,
                                                                       18517, 7852, 18847, 1954,
                                                                       2020, 8974, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 20365, 0, 3,
                                                                       19177, 8578, 19573, 2152,
                                                                       2230, 9640, ncols, gamma,
                                                                       p, q);

                    compute_prim_osd_three_center_electron_repulsion_0(buffer, 20833, 0, 3,
                                                                       19573, 8776, 19969, 2230,
                                                                       2308, 9874, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsd_three_center_electron_repulsion_0(buffer, 21301, 0, 3,
                                                                       20365, 9640, 20833, 2464,
                                                                       2555, 10654, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21847, 3, 2737,
                                                                       2740, 10927, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21857, 3, 2740,
                                                                       2743, 10933, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21867, 3, 2743,
                                                                       2746, 10939, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21877, 3, 2746,
                                                                       2749, 10945, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21887, 3, 2749,
                                                                       2752, 10951, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21897, 3, 2752,
                                                                       2755, 10957, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21907, 3, 2755,
                                                                       2758, 10963, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21917, 3, 2758,
                                                                       2761, 10969, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21927, 3, 2761,
                                                                       2764, 10975, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21937, 3, 2764,
                                                                       2767, 10981, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21947, 3, 2767,
                                                                       2770, 10987, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21957, 3, 2770,
                                                                       2773, 10993, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 21967, 3, 2773,
                                                                       2776, 10999, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 21977, 0, 3,
                                                                       21847, 10927, 21857, 2782,
                                                                       2791, 11005, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22007, 0, 3,
                                                                       21857, 10933, 21867, 2791,
                                                                       2800, 11023, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22037, 0, 3,
                                                                       21867, 10939, 21877, 2800,
                                                                       2809, 11041, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22067, 0, 3,
                                                                       21877, 10945, 21887, 2809,
                                                                       2818, 11059, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22097, 0, 3,
                                                                       21887, 10951, 21897, 2818,
                                                                       2827, 11077, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22127, 0, 3,
                                                                       21897, 10957, 21907, 2827,
                                                                       2836, 11095, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22157, 0, 3,
                                                                       21907, 10963, 21917, 2836,
                                                                       2845, 11113, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22187, 0, 3,
                                                                       21917, 10969, 21927, 2845,
                                                                       2854, 11131, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22217, 0, 3,
                                                                       21927, 10975, 21937, 2854,
                                                                       2863, 11149, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22247, 0, 3,
                                                                       21937, 10981, 21947, 2863,
                                                                       2872, 11167, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22277, 0, 3,
                                                                       21947, 10987, 21957, 2872,
                                                                       2881, 11185, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 22307, 0, 3,
                                                                       21957, 10993, 21967, 2881,
                                                                       2890, 11203, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22337, 0, 3,
                                                                       21977, 11005, 22007, 2908,
                                                                       2926, 11221, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22397, 0, 3,
                                                                       22007, 11023, 22037, 2926,
                                                                       2944, 11257, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22457, 0, 3,
                                                                       22037, 11041, 22067, 2944,
                                                                       2962, 11293, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22517, 0, 3,
                                                                       22067, 11059, 22097, 2962,
                                                                       2980, 11329, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22577, 0, 3,
                                                                       22097, 11077, 22127, 2980,
                                                                       2998, 11365, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22637, 0, 3,
                                                                       22127, 11095, 22157, 2998,
                                                                       3016, 11401, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22697, 0, 3,
                                                                       22157, 11113, 22187, 3016,
                                                                       3034, 11437, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22757, 0, 3,
                                                                       22187, 11131, 22217, 3034,
                                                                       3052, 11473, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22817, 0, 3,
                                                                       22217, 11149, 22247, 3052,
                                                                       3070, 11509, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22877, 0, 3,
                                                                       22247, 11167, 22277, 3070,
                                                                       3088, 11545, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 22937, 0, 3,
                                                                       22277, 11185, 22307, 3088,
                                                                       3106, 11581, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22997, 0, 3,
                                                                       22337, 11221, 22397, 3142,
                                                                       3172, 11617, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23097, 0, 3,
                                                                       22397, 11257, 22457, 3172,
                                                                       3202, 11677, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23197, 0, 3,
                                                                       22457, 11293, 22517, 3202,
                                                                       3232, 11737, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23297, 0, 3,
                                                                       22517, 11329, 22577, 3232,
                                                                       3262, 11797, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23397, 0, 3,
                                                                       22577, 11365, 22637, 3262,
                                                                       3292, 11857, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23497, 0, 3,
                                                                       22637, 11401, 22697, 3292,
                                                                       3322, 11917, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23597, 0, 3,
                                                                       22697, 11437, 22757, 3322,
                                                                       3352, 11977, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23697, 0, 3,
                                                                       22757, 11473, 22817, 3352,
                                                                       3382, 12037, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23797, 0, 3,
                                                                       22817, 11509, 22877, 3382,
                                                                       3412, 12097, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23897, 0, 3,
                                                                       22877, 11545, 22937, 3412,
                                                                       3442, 12157, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23997, 0, 3,
                                                                       22997, 11617, 23097, 3502,
                                                                       3547, 12217, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24147, 0, 3,
                                                                       23097, 11677, 23197, 3547,
                                                                       3592, 12307, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24297, 0, 3,
                                                                       23197, 11737, 23297, 3592,
                                                                       3637, 12397, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24447, 0, 3,
                                                                       23297, 11797, 23397, 3637,
                                                                       3682, 12487, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24597, 0, 3,
                                                                       23397, 11857, 23497, 3682,
                                                                       3727, 12577, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24747, 0, 3,
                                                                       23497, 11917, 23597, 3727,
                                                                       3772, 12667, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 24897, 0, 3,
                                                                       23597, 11977, 23697, 3772,
                                                                       3817, 12757, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 25047, 0, 3,
                                                                       23697, 12037, 23797, 3817,
                                                                       3862, 12847, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 25197, 0, 3,
                                                                       23797, 12097, 23897, 3862,
                                                                       3907, 12937, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25347, 0, 3,
                                                                       23997, 12217, 24147, 3997,
                                                                       4060, 13027, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25557, 0, 3,
                                                                       24147, 12307, 24297, 4060,
                                                                       4123, 13153, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25767, 0, 3,
                                                                       24297, 12397, 24447, 4123,
                                                                       4186, 13279, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 25977, 0, 3,
                                                                       24447, 12487, 24597, 4186,
                                                                       4249, 13405, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26187, 0, 3,
                                                                       24597, 12577, 24747, 4249,
                                                                       4312, 13531, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26397, 0, 3,
                                                                       24747, 12667, 24897, 4312,
                                                                       4375, 13657, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26607, 0, 3,
                                                                       24897, 12757, 25047, 4375,
                                                                       4438, 13783, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 26817, 0, 3,
                                                                       25047, 12847, 25197, 4438,
                                                                       4501, 13909, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27027, 0, 3,
                                                                       25347, 13027, 25557, 4627,
                                                                       4711, 14035, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27307, 0, 3,
                                                                       25557, 13153, 25767, 4711,
                                                                       4795, 14203, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27587, 0, 3,
                                                                       25767, 13279, 25977, 4795,
                                                                       4879, 14371, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 27867, 0, 3,
                                                                       25977, 13405, 26187, 4879,
                                                                       4963, 14539, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28147, 0, 3,
                                                                       26187, 13531, 26397, 4963,
                                                                       5047, 14707, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28427, 0, 3,
                                                                       26397, 13657, 26607, 5047,
                                                                       5131, 14875, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 28707, 0, 3,
                                                                       26607, 13783, 26817, 5131,
                                                                       5215, 15043, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 28987, 0, 3,
                                                                       27027, 14035, 27307, 5383,
                                                                       5491, 15211, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 29347, 0, 3,
                                                                       27307, 14203, 27587, 5491,
                                                                       5599, 15427, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 29707, 0, 3,
                                                                       27587, 14371, 27867, 5599,
                                                                       5707, 15643, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 30067, 0, 3,
                                                                       27867, 14539, 28147, 5707,
                                                                       5815, 15859, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 30427, 0, 3,
                                                                       28147, 14707, 28427, 5815,
                                                                       5923, 16075, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 30787, 0, 3,
                                                                       28427, 14875, 28707, 5923,
                                                                       6031, 16291, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 31147, 0, 3,
                                                                       28987, 15211, 29347, 6247,
                                                                       6382, 16507, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 31597, 0, 3,
                                                                       29347, 15427, 29707, 6382,
                                                                       6517, 16777, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 32047, 0, 3,
                                                                       29707, 15643, 30067, 6517,
                                                                       6652, 17047, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 32497, 0, 3,
                                                                       30067, 15859, 30427, 6652,
                                                                       6787, 17317, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 32947, 0, 3,
                                                                       30427, 16075, 30787, 6787,
                                                                       6922, 17587, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 33397, 0, 3,
                                                                       31147, 16507, 31597, 7192,
                                                                       7357, 17857, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 33947, 0, 3,
                                                                       31597, 16777, 32047, 7357,
                                                                       7522, 18187, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 34497, 0, 3,
                                                                       32047, 17047, 32497, 7522,
                                                                       7687, 18517, ncols, gamma,
                                                                       p, q);

                    compute_prim_msf_three_center_electron_repulsion_0(buffer, 35047, 0, 3,
                                                                       32497, 17317, 32947, 7687,
                                                                       7852, 18847, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 35597, 0, 3,
                                                                       33397, 17857, 33947, 8182,
                                                                       8380, 19177, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 36257, 0, 3,
                                                                       33947, 18187, 34497, 8380,
                                                                       8578, 19573, ncols, gamma,
                                                                       p, q);

                    compute_prim_nsf_three_center_electron_repulsion_0(buffer, 36917, 0, 3,
                                                                       34497, 18517, 35047, 8578,
                                                                       8776, 19969, ncols, gamma,
                                                                       p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 37577, 0, 3,
                                                                       35597, 19177, 36257, 9172,
                                                                       9406, 20365, ncols, gamma,
                                                                       p, q);

                    compute_prim_osf_three_center_electron_repulsion_0(buffer, 38357, 0, 3,
                                                                       36257, 19573, 36917, 9406,
                                                                       9640, 20833, ncols, gamma,
                                                                       p, q);

                    compute_prim_qsf_three_center_electron_repulsion_0(buffer, 39137, 0, 3,
                                                                       37577, 20365, 38357,
                                                                       10108, 10381, 21301,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 40047, 27027, 280, ncols);

                    simdfunc::contract_primitives(buffer, 40523, 28987, 360, ncols);

                    simdfunc::contract_primitives(buffer, 41135, 31147, 450, ncols);

                    simdfunc::contract_primitives(buffer, 41900, 33397, 550, ncols);

                    simdfunc::contract_primitives(buffer, 42835, 35597, 660, ncols);

                    simdfunc::contract_primitives(buffer, 43957, 37577, 780, ncols);

                    simdfunc::contract_primitives(buffer, 45283, 39137, 910, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 40327, 40047, 28, 1, nmax);

        simdtrf::transform_f_inner(buffer, 40883, 40523, 36, 1, nmax);

        simdtrf::transform_f_inner(buffer, 41585, 41135, 45, 1, nmax);

        simdtrf::transform_f_inner(buffer, 42450, 41900, 55, 1, nmax);

        simdtrf::transform_f_inner(buffer, 43495, 42835, 66, 1, nmax);

        simdtrf::transform_f_inner(buffer, 44737, 43957, 78, 1, nmax);

        simdtrf::transform_f_inner(buffer, 46193, 45283, 91, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 46830, 40327, 40883, 7, nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 47418, 40883, 41585, 7, nmax);

        simdtrf::compute_hrr_lp_out_of_first(buffer, coordinates, 48174, 41585, 42450, 7, nmax);

        simdtrf::compute_hrr_mp_out_of_first(buffer, coordinates, 49119, 42450, 43495, 7, nmax);

        simdtrf::compute_hrr_np_out_of_first(buffer, coordinates, 50274, 43495, 44737, 7, nmax);

        simdtrf::compute_hrr_op_out_of_first(buffer, coordinates, 51660, 44737, 46193, 7, nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 53298, 46830, 47418, 7, nmax);

        simdtrf::compute_hrr_kd_out_of_first(buffer, coordinates, 54474, 47418, 48174, 7, nmax);

        simdtrf::compute_hrr_ld_out_of_first(buffer, coordinates, 55986, 48174, 49119, 7, nmax);

        simdtrf::compute_hrr_md_out_of_first(buffer, coordinates, 57876, 49119, 50274, 7, nmax);

        simdtrf::compute_hrr_nd_out_of_first(buffer, coordinates, 60186, 50274, 51660, 7, nmax);

        simdtrf::compute_hrr_if_out_of_first(buffer, coordinates, 62958, 53298, 54474, 7, nmax);

        simdtrf::compute_hrr_kf_out_of_first(buffer, coordinates, 64918, 54474, 55986, 7, nmax);

        simdtrf::compute_hrr_lf_out_of_first(buffer, coordinates, 67438, 55986, 57876, 7, nmax);

        simdtrf::compute_hrr_mf_out_of_first(buffer, coordinates, 70588, 57876, 60186, 7, nmax);

        simdtrf::compute_hrr_ig_out_of_first(buffer, coordinates, 74438, 62958, 64918, 7, nmax);

        simdtrf::compute_hrr_kg_out_of_first(buffer, coordinates, 77378, 64918, 67438, 7, nmax);

        simdtrf::compute_hrr_lg_out_of_first(buffer, coordinates, 81158, 67438, 70588, 7, nmax);

        simdtrf::compute_hrr_ih_out_of_first(buffer, coordinates, 85883, 74438, 77378, 7, nmax);

        simdtrf::compute_hrr_kh_out_of_first(buffer, coordinates, 89999, 77378, 81158, 7, nmax);

        simdtrf::compute_hrr_ii_out_of_first(buffer, coordinates, 95291, 85883, 89999, 7, nmax);

        simdtrf::transform_i_inner(buffer, 100779, 95291, 28, 7, nmax);

        simdtrf::transform_i_outer(values + n * npairs, nvalues, buffer, 100779, 91, nmax);
    }

    for (size_t m = 0; m < 1183; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
