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


#include "SimdThreeCenterElectronRepulsionRecFFL.hpp"

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
#include "SimdTransferFD.hpp"
#include "SimdTransferFF.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGD.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_ffl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ffl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 68871, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 833 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 68871, 56497, 4112, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 889, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 892, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 895, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 898, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 901, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 904, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 907, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 910, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 913, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 916, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 919, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 922, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 925, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 928, 3, 9, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 937, 3, 10, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 946, 3, 11, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 955, 3, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 964, 3, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 973, 3, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 982, 3, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 991, 3, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1000, 3, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1009, 3, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1018, 3, 19, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1027, 3, 20, 61,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1036, 3, 28, 76,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1054, 3, 31, 82,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1072, 3, 34, 88,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1090, 3, 37, 94,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1108, 3, 40, 100,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1126, 3, 43, 106,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1144, 3, 46, 112,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1162, 3, 49, 118,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1180, 3, 52, 124,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1198, 3, 55, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1216, 3, 58, 136,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1234, 3, 76, 162,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1264, 3, 82, 172,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1294, 3, 88, 182,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1324, 3, 94, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1354, 3, 100, 202,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1384, 3, 106, 212,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1414, 3, 112, 222,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1444, 3, 118, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1474, 3, 124, 242,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1504, 3, 130, 252,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1534, 3, 162, 292,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1579, 3, 172, 307,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1624, 3, 182, 322,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1669, 3, 192, 337,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1714, 3, 202, 352,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1759, 3, 212, 367,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1804, 3, 222, 382,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1849, 3, 232, 397,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1894, 3, 242, 412,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1939, 3, 292, 469,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2002, 3, 307, 490,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2065, 3, 322, 511,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2128, 3, 337, 532,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2191, 3, 352, 553,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2254, 3, 367, 574,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2317, 3, 382, 595,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2380, 3, 397, 616,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2443, 3, 469, 693,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2527, 3, 490, 721,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2611, 3, 511, 749,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2695, 3, 532, 777,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2779, 3, 553, 805,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2863, 3, 574, 833,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 2947, 3, 595, 861,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3031, 3, 7, 8,
                                                                       889, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3037, 3, 8, 9,
                                                                       892, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3043, 3, 9, 10,
                                                                       895, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3049, 3, 10, 11,
                                                                       898, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3055, 3, 11, 12,
                                                                       901, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3061, 3, 12, 13,
                                                                       904, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3067, 3, 13, 14,
                                                                       907, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3073, 3, 14, 15,
                                                                       910, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3079, 3, 15, 16,
                                                                       913, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3085, 3, 16, 17,
                                                                       916, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3091, 3, 17, 18,
                                                                       919, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3097, 3, 18, 19,
                                                                       922, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3103, 3, 19, 20,
                                                                       925, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3109, 0, 3, 3031,
                                                                       889, 3037, 928, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3127, 0, 3, 3037,
                                                                       892, 3043, 937, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3145, 0, 3, 3043,
                                                                       895, 3049, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3163, 0, 3, 3049,
                                                                       898, 3055, 955, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3181, 0, 3, 3055,
                                                                       901, 3061, 964, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3199, 0, 3, 3061,
                                                                       904, 3067, 973, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3217, 0, 3, 3067,
                                                                       907, 3073, 982, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3235, 0, 3, 3073,
                                                                       910, 3079, 991, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3253, 0, 3, 3079,
                                                                       913, 3085, 1000, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3271, 0, 3, 3085,
                                                                       916, 3091, 1009, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3289, 0, 3, 3091,
                                                                       919, 3097, 1018, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3307, 0, 3, 3097,
                                                                       922, 3103, 1027, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3325, 0, 3, 3109,
                                                                       928, 3127, 64, 70, 1036,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3361, 0, 3, 3127,
                                                                       937, 3145, 70, 76, 1054,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3397, 0, 3, 3145,
                                                                       946, 3163, 76, 82, 1072,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3433, 0, 3, 3163,
                                                                       955, 3181, 82, 88, 1090,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3469, 0, 3, 3181,
                                                                       964, 3199, 88, 94, 1108,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3505, 0, 3, 3199,
                                                                       973, 3217, 94, 100, 1126,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3541, 0, 3, 3217,
                                                                       982, 3235, 100, 106, 1144,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3577, 0, 3, 3235,
                                                                       991, 3253, 106, 112, 1162,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3613, 0, 3, 3253,
                                                                       1000, 3271, 112, 118,
                                                                       1180, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3649, 0, 3, 3271,
                                                                       1009, 3289, 118, 124,
                                                                       1198, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3685, 0, 3, 3289,
                                                                       1018, 3307, 124, 130,
                                                                       1216, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3721, 0, 3, 3325,
                                                                       1036, 3361, 142, 152,
                                                                       1234, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3781, 0, 3, 3361,
                                                                       1054, 3397, 152, 162,
                                                                       1264, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3841, 0, 3, 3397,
                                                                       1072, 3433, 162, 172,
                                                                       1294, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3901, 0, 3, 3433,
                                                                       1090, 3469, 172, 182,
                                                                       1324, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3961, 0, 3, 3469,
                                                                       1108, 3505, 182, 192,
                                                                       1354, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4021, 0, 3, 3505,
                                                                       1126, 3541, 192, 202,
                                                                       1384, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4081, 0, 3, 3541,
                                                                       1144, 3577, 202, 212,
                                                                       1414, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4141, 0, 3, 3577,
                                                                       1162, 3613, 212, 222,
                                                                       1444, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4201, 0, 3, 3613,
                                                                       1180, 3649, 222, 232,
                                                                       1474, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4261, 0, 3, 3649,
                                                                       1198, 3685, 232, 242,
                                                                       1504, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4321, 0, 3, 3721,
                                                                       1234, 3781, 262, 277,
                                                                       1534, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4411, 0, 3, 3781,
                                                                       1264, 3841, 277, 292,
                                                                       1579, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4501, 0, 3, 3841,
                                                                       1294, 3901, 292, 307,
                                                                       1624, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4591, 0, 3, 3901,
                                                                       1324, 3961, 307, 322,
                                                                       1669, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4681, 0, 3, 3961,
                                                                       1354, 4021, 322, 337,
                                                                       1714, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4771, 0, 3, 4021,
                                                                       1384, 4081, 337, 352,
                                                                       1759, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4861, 0, 3, 4081,
                                                                       1414, 4141, 352, 367,
                                                                       1804, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4951, 0, 3, 4141,
                                                                       1444, 4201, 367, 382,
                                                                       1849, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5041, 0, 3, 4201,
                                                                       1474, 4261, 382, 397,
                                                                       1894, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5131, 0, 3, 4321,
                                                                       1534, 4411, 427, 448,
                                                                       1939, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5257, 0, 3, 4411,
                                                                       1579, 4501, 448, 469,
                                                                       2002, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5383, 0, 3, 4501,
                                                                       1624, 4591, 469, 490,
                                                                       2065, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5509, 0, 3, 4591,
                                                                       1669, 4681, 490, 511,
                                                                       2128, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5635, 0, 3, 4681,
                                                                       1714, 4771, 511, 532,
                                                                       2191, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5761, 0, 3, 4771,
                                                                       1759, 4861, 532, 553,
                                                                       2254, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 5887, 0, 3, 4861,
                                                                       1804, 4951, 553, 574,
                                                                       2317, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6013, 0, 3, 4951,
                                                                       1849, 5041, 574, 595,
                                                                       2380, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6139, 0, 3, 5131,
                                                                       1939, 5257, 637, 665,
                                                                       2443, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6307, 0, 3, 5257,
                                                                       2002, 5383, 665, 693,
                                                                       2527, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6475, 0, 3, 5383,
                                                                       2065, 5509, 693, 721,
                                                                       2611, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6643, 0, 3, 5509,
                                                                       2128, 5635, 721, 749,
                                                                       2695, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6811, 0, 3, 5635,
                                                                       2191, 5761, 749, 777,
                                                                       2779, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 6979, 0, 3, 5761,
                                                                       2254, 5887, 777, 805,
                                                                       2863, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7147, 0, 3, 5887,
                                                                       2317, 6013, 805, 833,
                                                                       2947, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7315, 3, 889, 892,
                                                                       3043, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7325, 3, 892, 895,
                                                                       3049, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7335, 3, 895, 898,
                                                                       3055, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7345, 3, 898, 901,
                                                                       3061, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7355, 3, 901, 904,
                                                                       3067, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7365, 3, 904, 907,
                                                                       3073, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7375, 3, 907, 910,
                                                                       3079, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7385, 3, 910, 913,
                                                                       3085, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7395, 3, 913, 916,
                                                                       3091, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7405, 3, 916, 919,
                                                                       3097, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7415, 3, 919, 922,
                                                                       3103, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7425, 0, 3, 7315,
                                                                       3043, 7325, 928, 937,
                                                                       3145, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7455, 0, 3, 7325,
                                                                       3049, 7335, 937, 946,
                                                                       3163, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7485, 0, 3, 7335,
                                                                       3055, 7345, 946, 955,
                                                                       3181, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7515, 0, 3, 7345,
                                                                       3061, 7355, 955, 964,
                                                                       3199, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7545, 0, 3, 7355,
                                                                       3067, 7365, 964, 973,
                                                                       3217, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7575, 0, 3, 7365,
                                                                       3073, 7375, 973, 982,
                                                                       3235, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7605, 0, 3, 7375,
                                                                       3079, 7385, 982, 991,
                                                                       3253, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7635, 0, 3, 7385,
                                                                       3085, 7395, 991, 1000,
                                                                       3271, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7665, 0, 3, 7395,
                                                                       3091, 7405, 1000, 1009,
                                                                       3289, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7695, 0, 3, 7405,
                                                                       3097, 7415, 1009, 1018,
                                                                       3307, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7725, 0, 3, 7425,
                                                                       3145, 7455, 1036, 1054,
                                                                       3397, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7785, 0, 3, 7455,
                                                                       3163, 7485, 1054, 1072,
                                                                       3433, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7845, 0, 3, 7485,
                                                                       3181, 7515, 1072, 1090,
                                                                       3469, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7905, 0, 3, 7515,
                                                                       3199, 7545, 1090, 1108,
                                                                       3505, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7965, 0, 3, 7545,
                                                                       3217, 7575, 1108, 1126,
                                                                       3541, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8025, 0, 3, 7575,
                                                                       3235, 7605, 1126, 1144,
                                                                       3577, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8085, 0, 3, 7605,
                                                                       3253, 7635, 1144, 1162,
                                                                       3613, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8145, 0, 3, 7635,
                                                                       3271, 7665, 1162, 1180,
                                                                       3649, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8205, 0, 3, 7665,
                                                                       3289, 7695, 1180, 1198,
                                                                       3685, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8265, 0, 3, 7725,
                                                                       3397, 7785, 1234, 1264,
                                                                       3841, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8365, 0, 3, 7785,
                                                                       3433, 7845, 1264, 1294,
                                                                       3901, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8465, 0, 3, 7845,
                                                                       3469, 7905, 1294, 1324,
                                                                       3961, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8565, 0, 3, 7905,
                                                                       3505, 7965, 1324, 1354,
                                                                       4021, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8665, 0, 3, 7965,
                                                                       3541, 8025, 1354, 1384,
                                                                       4081, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8765, 0, 3, 8025,
                                                                       3577, 8085, 1384, 1414,
                                                                       4141, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8865, 0, 3, 8085,
                                                                       3613, 8145, 1414, 1444,
                                                                       4201, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8965, 0, 3, 8145,
                                                                       3649, 8205, 1444, 1474,
                                                                       4261, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9065, 0, 3, 8265,
                                                                       3841, 8365, 1534, 1579,
                                                                       4501, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9215, 0, 3, 8365,
                                                                       3901, 8465, 1579, 1624,
                                                                       4591, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9365, 0, 3, 8465,
                                                                       3961, 8565, 1624, 1669,
                                                                       4681, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9515, 0, 3, 8565,
                                                                       4021, 8665, 1669, 1714,
                                                                       4771, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9665, 0, 3, 8665,
                                                                       4081, 8765, 1714, 1759,
                                                                       4861, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9815, 0, 3, 8765,
                                                                       4141, 8865, 1759, 1804,
                                                                       4951, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9965, 0, 3, 8865,
                                                                       4201, 8965, 1804, 1849,
                                                                       5041, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 10115, 0, 3, 9065,
                                                                       4501, 9215, 1939, 2002,
                                                                       5383, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 10325, 0, 3, 9215,
                                                                       4591, 9365, 2002, 2065,
                                                                       5509, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 10535, 0, 3, 9365,
                                                                       4681, 9515, 2065, 2128,
                                                                       5635, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 10745, 0, 3, 9515,
                                                                       4771, 9665, 2128, 2191,
                                                                       5761, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 10955, 0, 3, 9665,
                                                                       4861, 9815, 2191, 2254,
                                                                       5887, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 11165, 0, 3, 9815,
                                                                       4951, 9965, 2254, 2317,
                                                                       6013, ncols, gamma, p,
                                                                       q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 11375, 0, 3,
                                                                       10115, 5383, 10325, 2443,
                                                                       2527, 6475, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 11655, 0, 3,
                                                                       10325, 5509, 10535, 2527,
                                                                       2611, 6643, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 11935, 0, 3,
                                                                       10535, 5635, 10745, 2611,
                                                                       2695, 6811, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 12215, 0, 3,
                                                                       10745, 5761, 10955, 2695,
                                                                       2779, 6979, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 12495, 0, 3,
                                                                       10955, 5887, 11165, 2779,
                                                                       2863, 7147, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12775, 3, 3031,
                                                                       3037, 7315, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12790, 3, 3037,
                                                                       3043, 7325, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12805, 3, 3043,
                                                                       3049, 7335, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12820, 3, 3049,
                                                                       3055, 7345, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12835, 3, 3055,
                                                                       3061, 7355, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12850, 3, 3061,
                                                                       3067, 7365, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12865, 3, 3067,
                                                                       3073, 7375, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12880, 3, 3073,
                                                                       3079, 7385, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12895, 3, 3079,
                                                                       3085, 7395, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12910, 3, 3085,
                                                                       3091, 7405, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 12925, 3, 3091,
                                                                       3097, 7415, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 12940, 0, 3,
                                                                       12775, 7315, 12790, 3109,
                                                                       3127, 7425, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 12985, 0, 3,
                                                                       12790, 7325, 12805, 3127,
                                                                       3145, 7455, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 13030, 0, 3,
                                                                       12805, 7335, 12820, 3145,
                                                                       3163, 7485, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 13075, 0, 3,
                                                                       12820, 7345, 12835, 3163,
                                                                       3181, 7515, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 13120, 0, 3,
                                                                       12835, 7355, 12850, 3181,
                                                                       3199, 7545, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 13165, 0, 3,
                                                                       12850, 7365, 12865, 3199,
                                                                       3217, 7575, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 13210, 0, 3,
                                                                       12865, 7375, 12880, 3217,
                                                                       3235, 7605, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 13255, 0, 3,
                                                                       12880, 7385, 12895, 3235,
                                                                       3253, 7635, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 13300, 0, 3,
                                                                       12895, 7395, 12910, 3253,
                                                                       3271, 7665, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 13345, 0, 3,
                                                                       12910, 7405, 12925, 3271,
                                                                       3289, 7695, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 13390, 0, 3,
                                                                       12940, 7425, 12985, 3325,
                                                                       3361, 7725, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 13480, 0, 3,
                                                                       12985, 7455, 13030, 3361,
                                                                       3397, 7785, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 13570, 0, 3,
                                                                       13030, 7485, 13075, 3397,
                                                                       3433, 7845, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 13660, 0, 3,
                                                                       13075, 7515, 13120, 3433,
                                                                       3469, 7905, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 13750, 0, 3,
                                                                       13120, 7545, 13165, 3469,
                                                                       3505, 7965, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 13840, 0, 3,
                                                                       13165, 7575, 13210, 3505,
                                                                       3541, 8025, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 13930, 0, 3,
                                                                       13210, 7605, 13255, 3541,
                                                                       3577, 8085, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14020, 0, 3,
                                                                       13255, 7635, 13300, 3577,
                                                                       3613, 8145, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14110, 0, 3,
                                                                       13300, 7665, 13345, 3613,
                                                                       3649, 8205, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 14200, 0, 3,
                                                                       13390, 7725, 13480, 3721,
                                                                       3781, 8265, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 14350, 0, 3,
                                                                       13480, 7785, 13570, 3781,
                                                                       3841, 8365, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 14500, 0, 3,
                                                                       13570, 7845, 13660, 3841,
                                                                       3901, 8465, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 14650, 0, 3,
                                                                       13660, 7905, 13750, 3901,
                                                                       3961, 8565, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 14800, 0, 3,
                                                                       13750, 7965, 13840, 3961,
                                                                       4021, 8665, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 14950, 0, 3,
                                                                       13840, 8025, 13930, 4021,
                                                                       4081, 8765, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15100, 0, 3,
                                                                       13930, 8085, 14020, 4081,
                                                                       4141, 8865, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15250, 0, 3,
                                                                       14020, 8145, 14110, 4141,
                                                                       4201, 8965, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 15400, 0, 3,
                                                                       14200, 8265, 14350, 4321,
                                                                       4411, 9065, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 15625, 0, 3,
                                                                       14350, 8365, 14500, 4411,
                                                                       4501, 9215, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 15850, 0, 3,
                                                                       14500, 8465, 14650, 4501,
                                                                       4591, 9365, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 16075, 0, 3,
                                                                       14650, 8565, 14800, 4591,
                                                                       4681, 9515, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 16300, 0, 3,
                                                                       14800, 8665, 14950, 4681,
                                                                       4771, 9665, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 16525, 0, 3,
                                                                       14950, 8765, 15100, 4771,
                                                                       4861, 9815, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 16750, 0, 3,
                                                                       15100, 8865, 15250, 4861,
                                                                       4951, 9965, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 16975, 0, 3,
                                                                       15400, 9065, 15625, 5131,
                                                                       5257, 10115, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 17290, 0, 3,
                                                                       15625, 9215, 15850, 5257,
                                                                       5383, 10325, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 17605, 0, 3,
                                                                       15850, 9365, 16075, 5383,
                                                                       5509, 10535, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 17920, 0, 3,
                                                                       16075, 9515, 16300, 5509,
                                                                       5635, 10745, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 18235, 0, 3,
                                                                       16300, 9665, 16525, 5635,
                                                                       5761, 10955, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 18550, 0, 3,
                                                                       16525, 9815, 16750, 5761,
                                                                       5887, 11165, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 18865, 0, 3,
                                                                       16975, 10115, 17290, 6139,
                                                                       6307, 11375, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 19285, 0, 3,
                                                                       17290, 10325, 17605, 6307,
                                                                       6475, 11655, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 19705, 0, 3,
                                                                       17605, 10535, 17920, 6475,
                                                                       6643, 11935, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 20125, 0, 3,
                                                                       17920, 10745, 18235, 6643,
                                                                       6811, 12215, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 20545, 0, 3,
                                                                       18235, 10955, 18550, 6811,
                                                                       6979, 12495, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 20965, 3, 7315,
                                                                       7325, 12805, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 20986, 3, 7325,
                                                                       7335, 12820, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21007, 3, 7335,
                                                                       7345, 12835, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21028, 3, 7345,
                                                                       7355, 12850, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21049, 3, 7355,
                                                                       7365, 12865, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21070, 3, 7365,
                                                                       7375, 12880, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21091, 3, 7375,
                                                                       7385, 12895, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21112, 3, 7385,
                                                                       7395, 12910, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 21133, 3, 7395,
                                                                       7405, 12925, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21154, 0, 3,
                                                                       20965, 12805, 20986, 7425,
                                                                       7455, 13030, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21217, 0, 3,
                                                                       20986, 12820, 21007, 7455,
                                                                       7485, 13075, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21280, 0, 3,
                                                                       21007, 12835, 21028, 7485,
                                                                       7515, 13120, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21343, 0, 3,
                                                                       21028, 12850, 21049, 7515,
                                                                       7545, 13165, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21406, 0, 3,
                                                                       21049, 12865, 21070, 7545,
                                                                       7575, 13210, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21469, 0, 3,
                                                                       21070, 12880, 21091, 7575,
                                                                       7605, 13255, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21532, 0, 3,
                                                                       21091, 12895, 21112, 7605,
                                                                       7635, 13300, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21595, 0, 3,
                                                                       21112, 12910, 21133, 7635,
                                                                       7665, 13345, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 21658, 0, 3,
                                                                       21154, 13030, 21217, 7725,
                                                                       7785, 13570, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 21784, 0, 3,
                                                                       21217, 13075, 21280, 7785,
                                                                       7845, 13660, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 21910, 0, 3,
                                                                       21280, 13120, 21343, 7845,
                                                                       7905, 13750, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 22036, 0, 3,
                                                                       21343, 13165, 21406, 7905,
                                                                       7965, 13840, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 22162, 0, 3,
                                                                       21406, 13210, 21469, 7965,
                                                                       8025, 13930, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 22288, 0, 3,
                                                                       21469, 13255, 21532, 8025,
                                                                       8085, 14020, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 22414, 0, 3,
                                                                       21532, 13300, 21595, 8085,
                                                                       8145, 14110, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 22540, 0, 3,
                                                                       21658, 13570, 21784, 8265,
                                                                       8365, 14500, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 22750, 0, 3,
                                                                       21784, 13660, 21910, 8365,
                                                                       8465, 14650, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 22960, 0, 3,
                                                                       21910, 13750, 22036, 8465,
                                                                       8565, 14800, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 23170, 0, 3,
                                                                       22036, 13840, 22162, 8565,
                                                                       8665, 14950, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 23380, 0, 3,
                                                                       22162, 13930, 22288, 8665,
                                                                       8765, 15100, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 23590, 0, 3,
                                                                       22288, 14020, 22414, 8765,
                                                                       8865, 15250, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 23800, 0, 3,
                                                                       22540, 14500, 22750, 9065,
                                                                       9215, 15850, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 24115, 0, 3,
                                                                       22750, 14650, 22960, 9215,
                                                                       9365, 16075, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 24430, 0, 3,
                                                                       22960, 14800, 23170, 9365,
                                                                       9515, 16300, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 24745, 0, 3,
                                                                       23170, 14950, 23380, 9515,
                                                                       9665, 16525, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 25060, 0, 3,
                                                                       23380, 15100, 23590, 9665,
                                                                       9815, 16750, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 25375, 0, 3,
                                                                       23800, 15850, 24115,
                                                                       10115, 10325, 17605,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 25816, 0, 3,
                                                                       24115, 16075, 24430,
                                                                       10325, 10535, 17920,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 26257, 0, 3,
                                                                       24430, 16300, 24745,
                                                                       10535, 10745, 18235,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 26698, 0, 3,
                                                                       24745, 16525, 25060,
                                                                       10745, 10955, 18550,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 27139, 0, 3,
                                                                       25375, 17605, 25816,
                                                                       11375, 11655, 19705,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 27727, 0, 3,
                                                                       25816, 17920, 26257,
                                                                       11655, 11935, 20125,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 28315, 0, 3,
                                                                       26257, 18235, 26698,
                                                                       11935, 12215, 20545,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 28903, 3, 12775,
                                                                       12790, 20965, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 28931, 3, 12790,
                                                                       12805, 20986, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 28959, 3, 12805,
                                                                       12820, 21007, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 28987, 3, 12820,
                                                                       12835, 21028, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29015, 3, 12835,
                                                                       12850, 21049, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29043, 3, 12850,
                                                                       12865, 21070, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29071, 3, 12865,
                                                                       12880, 21091, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29099, 3, 12880,
                                                                       12895, 21112, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 29127, 3, 12895,
                                                                       12910, 21133, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29155, 0, 3,
                                                                       28903, 20965, 28931,
                                                                       12940, 12985, 21154,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29239, 0, 3,
                                                                       28931, 20986, 28959,
                                                                       12985, 13030, 21217,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29323, 0, 3,
                                                                       28959, 21007, 28987,
                                                                       13030, 13075, 21280,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29407, 0, 3,
                                                                       28987, 21028, 29015,
                                                                       13075, 13120, 21343,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29491, 0, 3,
                                                                       29015, 21049, 29043,
                                                                       13120, 13165, 21406,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29575, 0, 3,
                                                                       29043, 21070, 29071,
                                                                       13165, 13210, 21469,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29659, 0, 3,
                                                                       29071, 21091, 29099,
                                                                       13210, 13255, 21532,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 29743, 0, 3,
                                                                       29099, 21112, 29127,
                                                                       13255, 13300, 21595,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 29827, 0, 3,
                                                                       29155, 21154, 29239,
                                                                       13390, 13480, 21658,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 29995, 0, 3,
                                                                       29239, 21217, 29323,
                                                                       13480, 13570, 21784,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 30163, 0, 3,
                                                                       29323, 21280, 29407,
                                                                       13570, 13660, 21910,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 30331, 0, 3,
                                                                       29407, 21343, 29491,
                                                                       13660, 13750, 22036,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 30499, 0, 3,
                                                                       29491, 21406, 29575,
                                                                       13750, 13840, 22162,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 30667, 0, 3,
                                                                       29575, 21469, 29659,
                                                                       13840, 13930, 22288,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 30835, 0, 3,
                                                                       29659, 21532, 29743,
                                                                       13930, 14020, 22414,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 31003, 0, 3,
                                                                       29827, 21658, 29995,
                                                                       14200, 14350, 22540,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 31283, 0, 3,
                                                                       29995, 21784, 30163,
                                                                       14350, 14500, 22750,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 31563, 0, 3,
                                                                       30163, 21910, 30331,
                                                                       14500, 14650, 22960,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 31843, 0, 3,
                                                                       30331, 22036, 30499,
                                                                       14650, 14800, 23170,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 32123, 0, 3,
                                                                       30499, 22162, 30667,
                                                                       14800, 14950, 23380,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 32403, 0, 3,
                                                                       30667, 22288, 30835,
                                                                       14950, 15100, 23590,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 32683, 0, 3,
                                                                       31003, 22540, 31283,
                                                                       15400, 15625, 23800,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 33103, 0, 3,
                                                                       31283, 22750, 31563,
                                                                       15625, 15850, 24115,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 33523, 0, 3,
                                                                       31563, 22960, 31843,
                                                                       15850, 16075, 24430,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 33943, 0, 3,
                                                                       31843, 23170, 32123,
                                                                       16075, 16300, 24745,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 34363, 0, 3,
                                                                       32123, 23380, 32403,
                                                                       16300, 16525, 25060,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 34783, 0, 3,
                                                                       32683, 23800, 33103,
                                                                       16975, 17290, 25375,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 35371, 0, 3,
                                                                       33103, 24115, 33523,
                                                                       17290, 17605, 25816,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 35959, 0, 3,
                                                                       33523, 24430, 33943,
                                                                       17605, 17920, 26257,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 36547, 0, 3,
                                                                       33943, 24745, 34363,
                                                                       17920, 18235, 26698,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 37135, 0, 3,
                                                                       34783, 25375, 35371,
                                                                       18865, 19285, 27139,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 37919, 0, 3,
                                                                       35371, 25816, 35959,
                                                                       19285, 19705, 27727,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 38703, 0, 3,
                                                                       35959, 26257, 36547,
                                                                       19705, 20125, 28315,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39487, 3, 20965,
                                                                       20986, 28959, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39523, 3, 20986,
                                                                       21007, 28987, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39559, 3, 21007,
                                                                       21028, 29015, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39595, 3, 21028,
                                                                       21049, 29043, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39631, 3, 21049,
                                                                       21070, 29071, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39667, 3, 21070,
                                                                       21091, 29099, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 39703, 3, 21091,
                                                                       21112, 29127, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 39739, 0, 3,
                                                                       39487, 28959, 39523,
                                                                       21154, 21217, 29323,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 39847, 0, 3,
                                                                       39523, 28987, 39559,
                                                                       21217, 21280, 29407,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 39955, 0, 3,
                                                                       39559, 29015, 39595,
                                                                       21280, 21343, 29491,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 40063, 0, 3,
                                                                       39595, 29043, 39631,
                                                                       21343, 21406, 29575,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 40171, 0, 3,
                                                                       39631, 29071, 39667,
                                                                       21406, 21469, 29659,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 40279, 0, 3,
                                                                       39667, 29099, 39703,
                                                                       21469, 21532, 29743,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 40387, 0, 3,
                                                                       39739, 29323, 39847,
                                                                       21658, 21784, 30163,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 40603, 0, 3,
                                                                       39847, 29407, 39955,
                                                                       21784, 21910, 30331,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 40819, 0, 3,
                                                                       39955, 29491, 40063,
                                                                       21910, 22036, 30499,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 41035, 0, 3,
                                                                       40063, 29575, 40171,
                                                                       22036, 22162, 30667,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 41251, 0, 3,
                                                                       40171, 29659, 40279,
                                                                       22162, 22288, 30835,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 41467, 0, 3,
                                                                       40387, 30163, 40603,
                                                                       22540, 22750, 31563,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 41827, 0, 3,
                                                                       40603, 30331, 40819,
                                                                       22750, 22960, 31843,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 42187, 0, 3,
                                                                       40819, 30499, 41035,
                                                                       22960, 23170, 32123,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 42547, 0, 3,
                                                                       41035, 30667, 41251,
                                                                       23170, 23380, 32403,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 42907, 0, 3,
                                                                       41467, 31563, 41827,
                                                                       23800, 24115, 33523,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 43447, 0, 3,
                                                                       41827, 31843, 42187,
                                                                       24115, 24430, 33943,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 43987, 0, 3,
                                                                       42187, 32123, 42547,
                                                                       24430, 24745, 34363,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 44527, 0, 3,
                                                                       42907, 33523, 43447,
                                                                       25375, 25816, 35959,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsk_three_center_electron_repulsion_0(buffer, 45283, 0, 3,
                                                                       43447, 33943, 43987,
                                                                       25816, 26257, 36547,
                                                                       ncols, gamma, p, q);

                    compute_prim_isk_three_center_electron_repulsion_0(buffer, 46039, 0, 3,
                                                                       44527, 35959, 45283,
                                                                       27139, 27727, 38703,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47047, 3, 28903,
                                                                       28931, 39487, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47092, 3, 28931,
                                                                       28959, 39523, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47137, 3, 28959,
                                                                       28987, 39559, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47182, 3, 28987,
                                                                       29015, 39595, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47227, 3, 29015,
                                                                       29043, 39631, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47272, 3, 29043,
                                                                       29071, 39667, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 47317, 3, 29071,
                                                                       29099, 39703, ncols,
                                                                       gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 47362, 0, 3,
                                                                       47047, 39487, 47092,
                                                                       29155, 29239, 39739,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 47497, 0, 3,
                                                                       47092, 39523, 47137,
                                                                       29239, 29323, 39847,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 47632, 0, 3,
                                                                       47137, 39559, 47182,
                                                                       29323, 29407, 39955,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 47767, 0, 3,
                                                                       47182, 39595, 47227,
                                                                       29407, 29491, 40063,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 47902, 0, 3,
                                                                       47227, 39631, 47272,
                                                                       29491, 29575, 40171,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 48037, 0, 3,
                                                                       47272, 39667, 47317,
                                                                       29575, 29659, 40279,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 48172, 0, 3,
                                                                       47362, 39739, 47497,
                                                                       29827, 29995, 40387,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 48442, 0, 3,
                                                                       47497, 39847, 47632,
                                                                       29995, 30163, 40603,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 48712, 0, 3,
                                                                       47632, 39955, 47767,
                                                                       30163, 30331, 40819,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 48982, 0, 3,
                                                                       47767, 40063, 47902,
                                                                       30331, 30499, 41035,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 49252, 0, 3,
                                                                       47902, 40171, 48037,
                                                                       30499, 30667, 41251,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 49522, 0, 3,
                                                                       48172, 40387, 48442,
                                                                       31003, 31283, 41467,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 49972, 0, 3,
                                                                       48442, 40603, 48712,
                                                                       31283, 31563, 41827,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 50422, 0, 3,
                                                                       48712, 40819, 48982,
                                                                       31563, 31843, 42187,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 50872, 0, 3,
                                                                       48982, 41035, 49252,
                                                                       31843, 32123, 42547,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 51322, 0, 3,
                                                                       49522, 41467, 49972,
                                                                       32683, 33103, 42907,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 51997, 0, 3,
                                                                       49972, 41827, 50422,
                                                                       33103, 33523, 43447,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 52672, 0, 3,
                                                                       50422, 42187, 50872,
                                                                       33523, 33943, 43987,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 53347, 0, 3,
                                                                       51322, 42907, 51997,
                                                                       34783, 35371, 44527,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsl_three_center_electron_repulsion_0(buffer, 54292, 0, 3,
                                                                       51997, 43447, 52672,
                                                                       35371, 35959, 45283,
                                                                       ncols, gamma, p, q);

                    compute_prim_isl_three_center_electron_repulsion_0(buffer, 55237, 0, 3,
                                                                       53347, 44527, 54292,
                                                                       37135, 37919, 46039,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 56497, 49522, 450, ncols);

                    simdfunc::contract_primitives(buffer, 57117, 51322, 675, ncols);

                    simdfunc::contract_primitives(buffer, 58047, 53347, 945, ncols);

                    simdfunc::contract_primitives(buffer, 59349, 55237, 1260, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 56947, 56497, 10, 1, nmax);

        simdtrf::transform_l_inner(buffer, 57792, 57117, 15, 1, nmax);

        simdtrf::transform_l_inner(buffer, 58992, 58047, 21, 1, nmax);

        simdtrf::transform_l_inner(buffer, 60609, 59349, 28, 1, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 61085, 56947, 57792, 17,
                                             nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 61595, 57792, 58992, 17,
                                             nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 62360, 58992, 60609, 17,
                                             nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 63431, 61085, 61595, 17,
                                             nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 64451, 61595, 62360, 17,
                                             nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 65981, 63431, 64451, 17,
                                             nmax);

        simdtrf::transform_f_inner(buffer, 67681, 65981, 10, 17, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 67681, 119, nmax);
    }

    for (size_t m = 0; m < 833; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
