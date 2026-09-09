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


#include "SimdThreeCenterElectronRepulsionRecHFI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
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
#include "SimdTransferHD.hpp"
#include "SimdTransferHF.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferID.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransferKP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_hfi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_hfi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 73875, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1001 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 73875, 56767, 4745, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 14,
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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1492, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1495, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1498, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1501, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1504, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1507, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1510, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1513, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1516, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1519, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1522, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1525, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1528, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1531, 3, 9, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1540, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1549, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1558, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1567, 3, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1576, 3, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1585, 3, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1594, 3, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1603, 3, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1612, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1621, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1630, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1639, 3, 28, 76,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1657, 3, 31, 82,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1675, 3, 34, 88,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1693, 3, 37, 94,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1711, 3, 40, 100,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1729, 3, 43, 106,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1747, 3, 46, 112,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1765, 3, 49, 118,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1783, 3, 52, 124,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1801, 3, 55, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1819, 3, 58, 136,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1837, 3, 76, 162,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1867, 3, 82, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1897, 3, 88, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1927, 3, 94, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1957, 3, 100, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1987, 3, 106, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2017, 3, 112, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2047, 3, 118, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2077, 3, 124, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2107, 3, 130, 252,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2137, 3, 162, 292,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2182, 3, 172, 307,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2227, 3, 182, 322,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2272, 3, 192, 337,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2317, 3, 202, 352,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2362, 3, 212, 367,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2407, 3, 222, 382,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2452, 3, 232, 397,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2497, 3, 242, 412,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2542, 3, 292, 469,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2605, 3, 307, 490,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2668, 3, 322, 511,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2731, 3, 337, 532,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2794, 3, 352, 553,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2857, 3, 367, 574,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2920, 3, 382, 595,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2983, 3, 397, 616,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3046, 3, 469, 693,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3130, 3, 490, 721,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3214, 3, 511, 749,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3298, 3, 532, 777,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3382, 3, 553, 805,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3466, 3, 574, 833,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3550, 3, 595, 861,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3634, 3, 693, 961,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3742, 3, 721, 997,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3850, 3, 749,
                                                                       1033, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3958, 3, 777,
                                                                       1069, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4066, 3, 805,
                                                                       1105, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4174, 3, 833,
                                                                       1141, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4282, 3, 961,
                                                                       1267, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4417, 3, 997,
                                                                       1312, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4552, 3, 1033,
                                                                       1357, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4687, 3, 1069,
                                                                       1402, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4822, 3, 1105,
                                                                       1447, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4957, 3, 7, 8,
                                                                       1492, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4963, 3, 8, 9,
                                                                       1495, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4969, 3, 9, 10,
                                                                       1498, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4975, 3, 10, 11,
                                                                       1501, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4981, 3, 11, 12,
                                                                       1504, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4987, 3, 12, 13,
                                                                       1507, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4993, 3, 13, 14,
                                                                       1510, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4999, 3, 14, 15,
                                                                       1513, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5005, 3, 15, 16,
                                                                       1516, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5011, 3, 16, 17,
                                                                       1519, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5017, 3, 17, 18,
                                                                       1522, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5023, 3, 18, 19,
                                                                       1525, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5029, 3, 19, 20,
                                                                       1528, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5035, 0, 3, 4957,
                                                                       1492, 4963, 1531, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5053, 0, 3, 4963,
                                                                       1495, 4969, 1540, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5071, 0, 3, 4969,
                                                                       1498, 4975, 1549, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5089, 0, 3, 4975,
                                                                       1501, 4981, 1558, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5107, 0, 3, 4981,
                                                                       1504, 4987, 1567, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5125, 0, 3, 4987,
                                                                       1507, 4993, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5143, 0, 3, 4993,
                                                                       1510, 4999, 1585, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5161, 0, 3, 4999,
                                                                       1513, 5005, 1594, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5179, 0, 3, 5005,
                                                                       1516, 5011, 1603, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5197, 0, 3, 5011,
                                                                       1519, 5017, 1612, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5215, 0, 3, 5017,
                                                                       1522, 5023, 1621, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5233, 0, 3, 5023,
                                                                       1525, 5029, 1630, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5251, 0, 3, 5035,
                                                                       1531, 5053, 64, 70, 1639,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5287, 0, 3, 5053,
                                                                       1540, 5071, 70, 76, 1657,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5323, 0, 3, 5071,
                                                                       1549, 5089, 76, 82, 1675,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5359, 0, 3, 5089,
                                                                       1558, 5107, 82, 88, 1693,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5395, 0, 3, 5107,
                                                                       1567, 5125, 88, 94, 1711,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5431, 0, 3, 5125,
                                                                       1576, 5143, 94, 100, 1729,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5467, 0, 3, 5143,
                                                                       1585, 5161, 100, 106,
                                                                       1747, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5503, 0, 3, 5161,
                                                                       1594, 5179, 106, 112,
                                                                       1765, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5539, 0, 3, 5179,
                                                                       1603, 5197, 112, 118,
                                                                       1783, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5575, 0, 3, 5197,
                                                                       1612, 5215, 118, 124,
                                                                       1801, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5611, 0, 3, 5215,
                                                                       1621, 5233, 124, 130,
                                                                       1819, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5647, 0, 3, 5251,
                                                                       1639, 5287, 142, 152,
                                                                       1837, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5707, 0, 3, 5287,
                                                                       1657, 5323, 152, 162,
                                                                       1867, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5767, 0, 3, 5323,
                                                                       1675, 5359, 162, 172,
                                                                       1897, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5827, 0, 3, 5359,
                                                                       1693, 5395, 172, 182,
                                                                       1927, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5887, 0, 3, 5395,
                                                                       1711, 5431, 182, 192,
                                                                       1957, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5947, 0, 3, 5431,
                                                                       1729, 5467, 192, 202,
                                                                       1987, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6007, 0, 3, 5467,
                                                                       1747, 5503, 202, 212,
                                                                       2017, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6067, 0, 3, 5503,
                                                                       1765, 5539, 212, 222,
                                                                       2047, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6127, 0, 3, 5539,
                                                                       1783, 5575, 222, 232,
                                                                       2077, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6187, 0, 3, 5575,
                                                                       1801, 5611, 232, 242,
                                                                       2107, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6247, 0, 3, 5647,
                                                                       1837, 5707, 262, 277,
                                                                       2137, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6337, 0, 3, 5707,
                                                                       1867, 5767, 277, 292,
                                                                       2182, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6427, 0, 3, 5767,
                                                                       1897, 5827, 292, 307,
                                                                       2227, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6517, 0, 3, 5827,
                                                                       1927, 5887, 307, 322,
                                                                       2272, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6607, 0, 3, 5887,
                                                                       1957, 5947, 322, 337,
                                                                       2317, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6697, 0, 3, 5947,
                                                                       1987, 6007, 337, 352,
                                                                       2362, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6787, 0, 3, 6007,
                                                                       2017, 6067, 352, 367,
                                                                       2407, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6877, 0, 3, 6067,
                                                                       2047, 6127, 367, 382,
                                                                       2452, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6967, 0, 3, 6127,
                                                                       2077, 6187, 382, 397,
                                                                       2497, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7057, 0, 3, 6247,
                                                                       2137, 6337, 427, 448,
                                                                       2542, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7183, 0, 3, 6337,
                                                                       2182, 6427, 448, 469,
                                                                       2605, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7309, 0, 3, 6427,
                                                                       2227, 6517, 469, 490,
                                                                       2668, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7435, 0, 3, 6517,
                                                                       2272, 6607, 490, 511,
                                                                       2731, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7561, 0, 3, 6607,
                                                                       2317, 6697, 511, 532,
                                                                       2794, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7687, 0, 3, 6697,
                                                                       2362, 6787, 532, 553,
                                                                       2857, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7813, 0, 3, 6787,
                                                                       2407, 6877, 553, 574,
                                                                       2920, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7939, 0, 3, 6877,
                                                                       2452, 6967, 574, 595,
                                                                       2983, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8065, 0, 3, 7057,
                                                                       2542, 7183, 637, 665,
                                                                       3046, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8233, 0, 3, 7183,
                                                                       2605, 7309, 665, 693,
                                                                       3130, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8401, 0, 3, 7309,
                                                                       2668, 7435, 693, 721,
                                                                       3214, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8569, 0, 3, 7435,
                                                                       2731, 7561, 721, 749,
                                                                       3298, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8737, 0, 3, 7561,
                                                                       2794, 7687, 749, 777,
                                                                       3382, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8905, 0, 3, 7687,
                                                                       2857, 7813, 777, 805,
                                                                       3466, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 9073, 0, 3, 7813,
                                                                       2920, 7939, 805, 833,
                                                                       3550, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 9241, 0, 3, 8065,
                                                                       3046, 8233, 889, 925,
                                                                       3634, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 9457, 0, 3, 8233,
                                                                       3130, 8401, 925, 961,
                                                                       3742, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 9673, 0, 3, 8401,
                                                                       3214, 8569, 961, 997,
                                                                       3850, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 9889, 0, 3, 8569,
                                                                       3298, 8737, 997, 1033,
                                                                       3958, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 10105, 0, 3, 8737,
                                                                       3382, 8905, 1033, 1069,
                                                                       4066, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 10321, 0, 3, 8905,
                                                                       3466, 9073, 1069, 1105,
                                                                       4174, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 10537, 0, 3, 9241,
                                                                       3634, 9457, 1177, 1222,
                                                                       4282, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 10807, 0, 3, 9457,
                                                                       3742, 9673, 1222, 1267,
                                                                       4417, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 11077, 0, 3, 9673,
                                                                       3850, 9889, 1267, 1312,
                                                                       4552, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 11347, 0, 3, 9889,
                                                                       3958, 10105, 1312, 1357,
                                                                       4687, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 11617, 0, 3,
                                                                       10105, 4066, 10321, 1357,
                                                                       1402, 4822, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11887, 3, 1492,
                                                                       1495, 4969, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11897, 3, 1495,
                                                                       1498, 4975, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11907, 3, 1498,
                                                                       1501, 4981, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11917, 3, 1501,
                                                                       1504, 4987, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11927, 3, 1504,
                                                                       1507, 4993, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11937, 3, 1507,
                                                                       1510, 4999, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11947, 3, 1510,
                                                                       1513, 5005, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11957, 3, 1513,
                                                                       1516, 5011, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11967, 3, 1516,
                                                                       1519, 5017, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11977, 3, 1519,
                                                                       1522, 5023, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 11987, 3, 1522,
                                                                       1525, 5029, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 11997, 0, 3,
                                                                       11887, 4969, 11897, 1531,
                                                                       1540, 5071, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12027, 0, 3,
                                                                       11897, 4975, 11907, 1540,
                                                                       1549, 5089, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12057, 0, 3,
                                                                       11907, 4981, 11917, 1549,
                                                                       1558, 5107, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12087, 0, 3,
                                                                       11917, 4987, 11927, 1558,
                                                                       1567, 5125, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12117, 0, 3,
                                                                       11927, 4993, 11937, 1567,
                                                                       1576, 5143, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12147, 0, 3,
                                                                       11937, 4999, 11947, 1576,
                                                                       1585, 5161, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12177, 0, 3,
                                                                       11947, 5005, 11957, 1585,
                                                                       1594, 5179, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12207, 0, 3,
                                                                       11957, 5011, 11967, 1594,
                                                                       1603, 5197, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12237, 0, 3,
                                                                       11967, 5017, 11977, 1603,
                                                                       1612, 5215, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12267, 0, 3,
                                                                       11977, 5023, 11987, 1612,
                                                                       1621, 5233, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12297, 0, 3,
                                                                       11997, 5071, 12027, 1639,
                                                                       1657, 5323, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12357, 0, 3,
                                                                       12027, 5089, 12057, 1657,
                                                                       1675, 5359, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12417, 0, 3,
                                                                       12057, 5107, 12087, 1675,
                                                                       1693, 5395, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12477, 0, 3,
                                                                       12087, 5125, 12117, 1693,
                                                                       1711, 5431, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12537, 0, 3,
                                                                       12117, 5143, 12147, 1711,
                                                                       1729, 5467, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12597, 0, 3,
                                                                       12147, 5161, 12177, 1729,
                                                                       1747, 5503, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12657, 0, 3,
                                                                       12177, 5179, 12207, 1747,
                                                                       1765, 5539, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12717, 0, 3,
                                                                       12207, 5197, 12237, 1765,
                                                                       1783, 5575, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12777, 0, 3,
                                                                       12237, 5215, 12267, 1783,
                                                                       1801, 5611, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 12837, 0, 3,
                                                                       12297, 5323, 12357, 1837,
                                                                       1867, 5767, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 12937, 0, 3,
                                                                       12357, 5359, 12417, 1867,
                                                                       1897, 5827, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13037, 0, 3,
                                                                       12417, 5395, 12477, 1897,
                                                                       1927, 5887, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13137, 0, 3,
                                                                       12477, 5431, 12537, 1927,
                                                                       1957, 5947, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13237, 0, 3,
                                                                       12537, 5467, 12597, 1957,
                                                                       1987, 6007, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13337, 0, 3,
                                                                       12597, 5503, 12657, 1987,
                                                                       2017, 6067, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13437, 0, 3,
                                                                       12657, 5539, 12717, 2017,
                                                                       2047, 6127, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13537, 0, 3,
                                                                       12717, 5575, 12777, 2047,
                                                                       2077, 6187, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13637, 0, 3,
                                                                       12837, 5767, 12937, 2137,
                                                                       2182, 6427, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13787, 0, 3,
                                                                       12937, 5827, 13037, 2182,
                                                                       2227, 6517, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13937, 0, 3,
                                                                       13037, 5887, 13137, 2227,
                                                                       2272, 6607, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 14087, 0, 3,
                                                                       13137, 5947, 13237, 2272,
                                                                       2317, 6697, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 14237, 0, 3,
                                                                       13237, 6007, 13337, 2317,
                                                                       2362, 6787, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 14387, 0, 3,
                                                                       13337, 6067, 13437, 2362,
                                                                       2407, 6877, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 14537, 0, 3,
                                                                       13437, 6127, 13537, 2407,
                                                                       2452, 6967, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 14687, 0, 3,
                                                                       13637, 6427, 13787, 2542,
                                                                       2605, 7309, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 14897, 0, 3,
                                                                       13787, 6517, 13937, 2605,
                                                                       2668, 7435, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 15107, 0, 3,
                                                                       13937, 6607, 14087, 2668,
                                                                       2731, 7561, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 15317, 0, 3,
                                                                       14087, 6697, 14237, 2731,
                                                                       2794, 7687, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 15527, 0, 3,
                                                                       14237, 6787, 14387, 2794,
                                                                       2857, 7813, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 15737, 0, 3,
                                                                       14387, 6877, 14537, 2857,
                                                                       2920, 7939, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 15947, 0, 3,
                                                                       14687, 7309, 14897, 3046,
                                                                       3130, 8401, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 16227, 0, 3,
                                                                       14897, 7435, 15107, 3130,
                                                                       3214, 8569, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 16507, 0, 3,
                                                                       15107, 7561, 15317, 3214,
                                                                       3298, 8737, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 16787, 0, 3,
                                                                       15317, 7687, 15527, 3298,
                                                                       3382, 8905, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 17067, 0, 3,
                                                                       15527, 7813, 15737, 3382,
                                                                       3466, 9073, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 17347, 0, 3,
                                                                       15947, 8401, 16227, 3634,
                                                                       3742, 9673, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 17707, 0, 3,
                                                                       16227, 8569, 16507, 3742,
                                                                       3850, 9889, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 18067, 0, 3,
                                                                       16507, 8737, 16787, 3850,
                                                                       3958, 10105, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 18427, 0, 3,
                                                                       16787, 8905, 17067, 3958,
                                                                       4066, 10321, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 18787, 0, 3,
                                                                       17347, 9673, 17707, 4282,
                                                                       4417, 11077, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 19237, 0, 3,
                                                                       17707, 9889, 18067, 4417,
                                                                       4552, 11347, ncols, gamma,
                                                                       p, q);

                    compute_prim_lsf_three_center_electron_repulsion_0(buffer, 19687, 0, 3,
                                                                       18067, 10105, 18427, 4552,
                                                                       4687, 11617, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20137, 3, 4957,
                                                                       4963, 11887, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20152, 3, 4963,
                                                                       4969, 11897, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20167, 3, 4969,
                                                                       4975, 11907, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20182, 3, 4975,
                                                                       4981, 11917, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20197, 3, 4981,
                                                                       4987, 11927, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20212, 3, 4987,
                                                                       4993, 11937, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20227, 3, 4993,
                                                                       4999, 11947, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20242, 3, 4999,
                                                                       5005, 11957, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20257, 3, 5005,
                                                                       5011, 11967, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20272, 3, 5011,
                                                                       5017, 11977, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 20287, 3, 5017,
                                                                       5023, 11987, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20302, 0, 3,
                                                                       20137, 11887, 20152, 5035,
                                                                       5053, 11997, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20347, 0, 3,
                                                                       20152, 11897, 20167, 5053,
                                                                       5071, 12027, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20392, 0, 3,
                                                                       20167, 11907, 20182, 5071,
                                                                       5089, 12057, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20437, 0, 3,
                                                                       20182, 11917, 20197, 5089,
                                                                       5107, 12087, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20482, 0, 3,
                                                                       20197, 11927, 20212, 5107,
                                                                       5125, 12117, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20527, 0, 3,
                                                                       20212, 11937, 20227, 5125,
                                                                       5143, 12147, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20572, 0, 3,
                                                                       20227, 11947, 20242, 5143,
                                                                       5161, 12177, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20617, 0, 3,
                                                                       20242, 11957, 20257, 5161,
                                                                       5179, 12207, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20662, 0, 3,
                                                                       20257, 11967, 20272, 5179,
                                                                       5197, 12237, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 20707, 0, 3,
                                                                       20272, 11977, 20287, 5197,
                                                                       5215, 12267, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 20752, 0, 3,
                                                                       20302, 11997, 20347, 5251,
                                                                       5287, 12297, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 20842, 0, 3,
                                                                       20347, 12027, 20392, 5287,
                                                                       5323, 12357, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 20932, 0, 3,
                                                                       20392, 12057, 20437, 5323,
                                                                       5359, 12417, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21022, 0, 3,
                                                                       20437, 12087, 20482, 5359,
                                                                       5395, 12477, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21112, 0, 3,
                                                                       20482, 12117, 20527, 5395,
                                                                       5431, 12537, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21202, 0, 3,
                                                                       20527, 12147, 20572, 5431,
                                                                       5467, 12597, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21292, 0, 3,
                                                                       20572, 12177, 20617, 5467,
                                                                       5503, 12657, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21382, 0, 3,
                                                                       20617, 12207, 20662, 5503,
                                                                       5539, 12717, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21472, 0, 3,
                                                                       20662, 12237, 20707, 5539,
                                                                       5575, 12777, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 21562, 0, 3,
                                                                       20752, 12297, 20842, 5647,
                                                                       5707, 12837, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 21712, 0, 3,
                                                                       20842, 12357, 20932, 5707,
                                                                       5767, 12937, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 21862, 0, 3,
                                                                       20932, 12417, 21022, 5767,
                                                                       5827, 13037, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22012, 0, 3,
                                                                       21022, 12477, 21112, 5827,
                                                                       5887, 13137, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22162, 0, 3,
                                                                       21112, 12537, 21202, 5887,
                                                                       5947, 13237, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22312, 0, 3,
                                                                       21202, 12597, 21292, 5947,
                                                                       6007, 13337, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22462, 0, 3,
                                                                       21292, 12657, 21382, 6007,
                                                                       6067, 13437, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22612, 0, 3,
                                                                       21382, 12717, 21472, 6067,
                                                                       6127, 13537, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 22762, 0, 3,
                                                                       21562, 12837, 21712, 6247,
                                                                       6337, 13637, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 22987, 0, 3,
                                                                       21712, 12937, 21862, 6337,
                                                                       6427, 13787, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 23212, 0, 3,
                                                                       21862, 13037, 22012, 6427,
                                                                       6517, 13937, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 23437, 0, 3,
                                                                       22012, 13137, 22162, 6517,
                                                                       6607, 14087, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 23662, 0, 3,
                                                                       22162, 13237, 22312, 6607,
                                                                       6697, 14237, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 23887, 0, 3,
                                                                       22312, 13337, 22462, 6697,
                                                                       6787, 14387, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 24112, 0, 3,
                                                                       22462, 13437, 22612, 6787,
                                                                       6877, 14537, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 24337, 0, 3,
                                                                       22762, 13637, 22987, 7057,
                                                                       7183, 14687, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 24652, 0, 3,
                                                                       22987, 13787, 23212, 7183,
                                                                       7309, 14897, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 24967, 0, 3,
                                                                       23212, 13937, 23437, 7309,
                                                                       7435, 15107, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 25282, 0, 3,
                                                                       23437, 14087, 23662, 7435,
                                                                       7561, 15317, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 25597, 0, 3,
                                                                       23662, 14237, 23887, 7561,
                                                                       7687, 15527, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 25912, 0, 3,
                                                                       23887, 14387, 24112, 7687,
                                                                       7813, 15737, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 26227, 0, 3,
                                                                       24337, 14687, 24652, 8065,
                                                                       8233, 15947, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 26647, 0, 3,
                                                                       24652, 14897, 24967, 8233,
                                                                       8401, 16227, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 27067, 0, 3,
                                                                       24967, 15107, 25282, 8401,
                                                                       8569, 16507, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 27487, 0, 3,
                                                                       25282, 15317, 25597, 8569,
                                                                       8737, 16787, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 27907, 0, 3,
                                                                       25597, 15527, 25912, 8737,
                                                                       8905, 17067, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 28327, 0, 3,
                                                                       26227, 15947, 26647, 9241,
                                                                       9457, 17347, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 28867, 0, 3,
                                                                       26647, 16227, 27067, 9457,
                                                                       9673, 17707, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 29407, 0, 3,
                                                                       27067, 16507, 27487, 9673,
                                                                       9889, 18067, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 29947, 0, 3,
                                                                       27487, 16787, 27907, 9889,
                                                                       10105, 18427, ncols,
                                                                       gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 30487, 0, 3,
                                                                       28327, 17347, 28867,
                                                                       10537, 10807, 18787,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 31162, 0, 3,
                                                                       28867, 17707, 29407,
                                                                       10807, 11077, 19237,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsg_three_center_electron_repulsion_0(buffer, 31837, 0, 3,
                                                                       29407, 18067, 29947,
                                                                       11077, 11347, 19687,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32512, 3, 11887,
                                                                       11897, 20167, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32533, 3, 11897,
                                                                       11907, 20182, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32554, 3, 11907,
                                                                       11917, 20197, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32575, 3, 11917,
                                                                       11927, 20212, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32596, 3, 11927,
                                                                       11937, 20227, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32617, 3, 11937,
                                                                       11947, 20242, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32638, 3, 11947,
                                                                       11957, 20257, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32659, 3, 11957,
                                                                       11967, 20272, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 32680, 3, 11967,
                                                                       11977, 20287, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 32701, 0, 3,
                                                                       32512, 20167, 32533,
                                                                       11997, 12027, 20392,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 32764, 0, 3,
                                                                       32533, 20182, 32554,
                                                                       12027, 12057, 20437,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 32827, 0, 3,
                                                                       32554, 20197, 32575,
                                                                       12057, 12087, 20482,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 32890, 0, 3,
                                                                       32575, 20212, 32596,
                                                                       12087, 12117, 20527,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 32953, 0, 3,
                                                                       32596, 20227, 32617,
                                                                       12117, 12147, 20572,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 33016, 0, 3,
                                                                       32617, 20242, 32638,
                                                                       12147, 12177, 20617,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 33079, 0, 3,
                                                                       32638, 20257, 32659,
                                                                       12177, 12207, 20662,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 33142, 0, 3,
                                                                       32659, 20272, 32680,
                                                                       12207, 12237, 20707,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33205, 0, 3,
                                                                       32701, 20392, 32764,
                                                                       12297, 12357, 20932,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33331, 0, 3,
                                                                       32764, 20437, 32827,
                                                                       12357, 12417, 21022,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33457, 0, 3,
                                                                       32827, 20482, 32890,
                                                                       12417, 12477, 21112,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33583, 0, 3,
                                                                       32890, 20527, 32953,
                                                                       12477, 12537, 21202,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33709, 0, 3,
                                                                       32953, 20572, 33016,
                                                                       12537, 12597, 21292,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33835, 0, 3,
                                                                       33016, 20617, 33079,
                                                                       12597, 12657, 21382,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 33961, 0, 3,
                                                                       33079, 20662, 33142,
                                                                       12657, 12717, 21472,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 34087, 0, 3,
                                                                       33205, 20932, 33331,
                                                                       12837, 12937, 21862,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 34297, 0, 3,
                                                                       33331, 21022, 33457,
                                                                       12937, 13037, 22012,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 34507, 0, 3,
                                                                       33457, 21112, 33583,
                                                                       13037, 13137, 22162,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 34717, 0, 3,
                                                                       33583, 21202, 33709,
                                                                       13137, 13237, 22312,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 34927, 0, 3,
                                                                       33709, 21292, 33835,
                                                                       13237, 13337, 22462,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 35137, 0, 3,
                                                                       33835, 21382, 33961,
                                                                       13337, 13437, 22612,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 35347, 0, 3,
                                                                       34087, 21862, 34297,
                                                                       13637, 13787, 23212,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 35662, 0, 3,
                                                                       34297, 22012, 34507,
                                                                       13787, 13937, 23437,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 35977, 0, 3,
                                                                       34507, 22162, 34717,
                                                                       13937, 14087, 23662,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 36292, 0, 3,
                                                                       34717, 22312, 34927,
                                                                       14087, 14237, 23887,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 36607, 0, 3,
                                                                       34927, 22462, 35137,
                                                                       14237, 14387, 24112,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 36922, 0, 3,
                                                                       35347, 23212, 35662,
                                                                       14687, 14897, 24967,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 37363, 0, 3,
                                                                       35662, 23437, 35977,
                                                                       14897, 15107, 25282,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 37804, 0, 3,
                                                                       35977, 23662, 36292,
                                                                       15107, 15317, 25597,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 38245, 0, 3,
                                                                       36292, 23887, 36607,
                                                                       15317, 15527, 25912,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 38686, 0, 3,
                                                                       36922, 24967, 37363,
                                                                       15947, 16227, 27067,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 39274, 0, 3,
                                                                       37363, 25282, 37804,
                                                                       16227, 16507, 27487,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 39862, 0, 3,
                                                                       37804, 25597, 38245,
                                                                       16507, 16787, 27907,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 40450, 0, 3,
                                                                       38686, 27067, 39274,
                                                                       17347, 17707, 29407,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 41206, 0, 3,
                                                                       39274, 27487, 39862,
                                                                       17707, 18067, 29947,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsh_three_center_electron_repulsion_0(buffer, 41962, 0, 3,
                                                                       40450, 29407, 41206,
                                                                       18787, 19237, 31837,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 42907, 3, 20137,
                                                                       20152, 32512, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 42935, 3, 20152,
                                                                       20167, 32533, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 42963, 3, 20167,
                                                                       20182, 32554, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 42991, 3, 20182,
                                                                       20197, 32575, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43019, 3, 20197,
                                                                       20212, 32596, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43047, 3, 20212,
                                                                       20227, 32617, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43075, 3, 20227,
                                                                       20242, 32638, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43103, 3, 20242,
                                                                       20257, 32659, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 43131, 3, 20257,
                                                                       20272, 32680, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43159, 0, 3,
                                                                       42907, 32512, 42935,
                                                                       20302, 20347, 32701,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43243, 0, 3,
                                                                       42935, 32533, 42963,
                                                                       20347, 20392, 32764,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43327, 0, 3,
                                                                       42963, 32554, 42991,
                                                                       20392, 20437, 32827,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43411, 0, 3,
                                                                       42991, 32575, 43019,
                                                                       20437, 20482, 32890,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43495, 0, 3,
                                                                       43019, 32596, 43047,
                                                                       20482, 20527, 32953,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43579, 0, 3,
                                                                       43047, 32617, 43075,
                                                                       20527, 20572, 33016,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43663, 0, 3,
                                                                       43075, 32638, 43103,
                                                                       20572, 20617, 33079,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 43747, 0, 3,
                                                                       43103, 32659, 43131,
                                                                       20617, 20662, 33142,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 43831, 0, 3,
                                                                       43159, 32701, 43243,
                                                                       20752, 20842, 33205,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 43999, 0, 3,
                                                                       43243, 32764, 43327,
                                                                       20842, 20932, 33331,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 44167, 0, 3,
                                                                       43327, 32827, 43411,
                                                                       20932, 21022, 33457,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 44335, 0, 3,
                                                                       43411, 32890, 43495,
                                                                       21022, 21112, 33583,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 44503, 0, 3,
                                                                       43495, 32953, 43579,
                                                                       21112, 21202, 33709,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 44671, 0, 3,
                                                                       43579, 33016, 43663,
                                                                       21202, 21292, 33835,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 44839, 0, 3,
                                                                       43663, 33079, 43747,
                                                                       21292, 21382, 33961,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 45007, 0, 3,
                                                                       43831, 33205, 43999,
                                                                       21562, 21712, 34087,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 45287, 0, 3,
                                                                       43999, 33331, 44167,
                                                                       21712, 21862, 34297,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 45567, 0, 3,
                                                                       44167, 33457, 44335,
                                                                       21862, 22012, 34507,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 45847, 0, 3,
                                                                       44335, 33583, 44503,
                                                                       22012, 22162, 34717,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 46127, 0, 3,
                                                                       44503, 33709, 44671,
                                                                       22162, 22312, 34927,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 46407, 0, 3,
                                                                       44671, 33835, 44839,
                                                                       22312, 22462, 35137,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 46687, 0, 3,
                                                                       45007, 34087, 45287,
                                                                       22762, 22987, 35347,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 47107, 0, 3,
                                                                       45287, 34297, 45567,
                                                                       22987, 23212, 35662,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 47527, 0, 3,
                                                                       45567, 34507, 45847,
                                                                       23212, 23437, 35977,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 47947, 0, 3,
                                                                       45847, 34717, 46127,
                                                                       23437, 23662, 36292,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 48367, 0, 3,
                                                                       46127, 34927, 46407,
                                                                       23662, 23887, 36607,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 48787, 0, 3,
                                                                       46687, 35347, 47107,
                                                                       24337, 24652, 36922,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 49375, 0, 3,
                                                                       47107, 35662, 47527,
                                                                       24652, 24967, 37363,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 49963, 0, 3,
                                                                       47527, 35977, 47947,
                                                                       24967, 25282, 37804,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 50551, 0, 3,
                                                                       47947, 36292, 48367,
                                                                       25282, 25597, 38245,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 51139, 0, 3,
                                                                       48787, 36922, 49375,
                                                                       26227, 26647, 38686,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 51923, 0, 3,
                                                                       49375, 37363, 49963,
                                                                       26647, 27067, 39274,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 52707, 0, 3,
                                                                       49963, 37804, 50551,
                                                                       27067, 27487, 39862,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 53491, 0, 3,
                                                                       51139, 38686, 51923,
                                                                       28327, 28867, 40450,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 54499, 0, 3,
                                                                       51923, 39274, 52707,
                                                                       28867, 29407, 41206,
                                                                       ncols, gamma, p, q);

                    compute_prim_lsi_three_center_electron_repulsion_0(buffer, 55507, 0, 3,
                                                                       53491, 40450, 54499,
                                                                       30487, 31162, 41962,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 56767, 48787, 588, ncols);

                    simdfunc::contract_primitives(buffer, 57628, 51139, 784, ncols);

                    simdfunc::contract_primitives(buffer, 58776, 53491, 1008, ncols);

                    simdfunc::contract_primitives(buffer, 60252, 55507, 1260, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 57355, 56767, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 58412, 57628, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 59784, 58776, 36, 1, nmax);

        simdtrf::transform_i_inner(buffer, 61512, 60252, 45, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 62097, 57355, 58412, 13,
                                             nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 62916, 58412, 59784, 13,
                                             nmax);

        simdtrf::compute_hrr_kp_out_of_first(buffer, coordinates, 64008, 59784, 61512, 13,
                                             nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 65412, 62097, 62916, 13,
                                             nmax);

        simdtrf::compute_hrr_id_out_of_first(buffer, coordinates, 67050, 62916, 64008, 13,
                                             nmax);

        simdtrf::compute_hrr_hf_out_of_first(buffer, coordinates, 69234, 65412, 67050, 13,
                                             nmax);

        simdtrf::transform_f_inner(buffer, 71964, 69234, 21, 13, nmax);

        simdtrf::transform_h_outer(values + n * npairs, nvalues, buffer, 71964, 91, nmax);
    }

    for (size_t m = 0; m < 1001; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
