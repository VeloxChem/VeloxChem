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


#include "SimdThreeCenterElectronRepulsionRsRecDDL.hpp"

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
#include "SimdTransferDD.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_ddl_three_center_electron_repulsion(double               *values,
                                               const size_t          npairs,
                                               const size_t          natoms,
                                               const CBasisFunction &a_function,
                                               const CBasisFunction &b_function,
                                               const CBasisFunction &c_function,
                                               const CSimdMatrix    &coordinates,
                                               const CSimdMatrix    &c_coordinates,
                                               CSimdMatrix          &buffer,
                                               const double          omega,
                                               const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_ddl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 48518, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 850 * natoms * npairs, 0.0);

        return;
    }

    const auto pi = mathconst::pi_value();

    // NOTE: a row of the values spans every atom pair of every atom on the ket
    // side, so a kernel handed the block of one atom steps by this to reach the
    // next component -- which is what lets it be the kernel a two-center form
    // uses, unchanged.

    const auto nvalues = natoms * npairs;

    simdfunc::compute_pair_exponents(a_function, b_function, coordinates, nmax);

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 48518, 41308, 3589, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 12,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 20, 3, 12,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 76, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 79, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 85, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 88, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 91, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 94, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 97, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 100, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 106, 0, 3, 7, 8,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 8, 9,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 9, 10,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 124, 0, 3, 10, 11,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 130, 0, 3, 11, 12,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 136, 0, 3, 12, 13,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 13, 14,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 14, 15,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 154, 0, 3, 15, 16,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 160, 0, 3, 16, 17,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 166, 0, 3, 17, 18,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 21, 22,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 178, 0, 3, 22, 23,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 184, 0, 3, 23, 24,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 190, 0, 3, 24, 25,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 196, 0, 3, 25, 26,
                                                                       82, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 26, 27,
                                                                       85, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 208, 0, 3, 27, 28,
                                                                       88, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 214, 0, 3, 28, 29,
                                                                       91, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 220, 0, 3, 29, 30,
                                                                       94, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 226, 0, 3, 30, 31,
                                                                       97, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 232, 0, 3, 31, 32,
                                                                       100, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 34, 37,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 37, 40,
                                                                       112, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 40, 43,
                                                                       118, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 43, 46,
                                                                       124, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 46, 49,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 49, 52,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 52, 55,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 55, 58,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 58, 61,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 61, 64,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 70, 73,
                                                                       172, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 73, 76,
                                                                       178, 184, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 76, 79,
                                                                       184, 190, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 79, 82,
                                                                       190, 196, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 82, 85,
                                                                       196, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 85, 88,
                                                                       202, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 88, 91,
                                                                       208, 214, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 91, 94,
                                                                       214, 220, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 94, 97,
                                                                       220, 226, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 97,
                                                                       100, 226, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 106,
                                                                       112, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 453, 0, 3, 112,
                                                                       118, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 118,
                                                                       124, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 483, 0, 3, 124,
                                                                       130, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 498, 0, 3, 130,
                                                                       136, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 513, 0, 3, 136,
                                                                       142, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 528, 0, 3, 142,
                                                                       148, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 543, 0, 3, 148,
                                                                       154, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 558, 0, 3, 154,
                                                                       160, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 573, 0, 3, 172,
                                                                       178, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 588, 0, 3, 178,
                                                                       184, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 603, 0, 3, 184,
                                                                       190, 358, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 618, 0, 3, 190,
                                                                       196, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 633, 0, 3, 196,
                                                                       202, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 648, 0, 3, 202,
                                                                       208, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 663, 0, 3, 208,
                                                                       214, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 678, 0, 3, 214,
                                                                       220, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 693, 0, 3, 220,
                                                                       226, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 708, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 711, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 714, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 717, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 720, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 723, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 726, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 729, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 732, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 735, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 738, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 741, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 744, 3, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 747, 3, 25, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 750, 3, 26, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 753, 3, 27, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 756, 3, 28, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 759, 3, 29, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 762, 3, 30, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 765, 3, 31, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 768, 3, 32, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 771, 3, 33, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 774, 3, 9, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 783, 3, 10, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 792, 3, 11, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 801, 3, 12, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 810, 3, 13, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 819, 3, 14, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 828, 3, 15, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 837, 3, 16, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 846, 3, 17, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 855, 3, 18, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 864, 3, 23, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 873, 3, 24, 79,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 882, 3, 25, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 891, 3, 26, 85,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 900, 3, 27, 88,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 909, 3, 28, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 918, 3, 29, 94,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 927, 3, 30, 97,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 936, 3, 31, 100,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 945, 3, 32, 103,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 954, 3, 40, 118,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 972, 3, 43, 124,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 990, 3, 46, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1008, 3, 49, 136,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1026, 3, 52, 142,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1044, 3, 55, 148,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1062, 3, 58, 154,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1080, 3, 61, 160,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1098, 3, 64, 166,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1116, 3, 76, 184,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1134, 3, 79, 190,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1152, 3, 82, 196,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1170, 3, 85, 202,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1188, 3, 88, 208,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1206, 3, 91, 214,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1224, 3, 94, 220,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1242, 3, 97, 226,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1260, 3, 100, 232,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1278, 3, 118, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1308, 3, 124, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1338, 3, 130, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1368, 3, 136, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1398, 3, 142, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1428, 3, 148, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1458, 3, 154, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1488, 3, 160, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1518, 3, 184, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1548, 3, 190, 368,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1578, 3, 196, 378,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1608, 3, 202, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1638, 3, 208, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1668, 3, 214, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1698, 3, 220, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1728, 3, 226, 428,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1758, 3, 258, 468,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1803, 3, 268, 483,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1848, 3, 278, 498,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1893, 3, 288, 513,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1938, 3, 298, 528,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1983, 3, 308, 543,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2028, 3, 318, 558,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2073, 3, 358, 603,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2118, 3, 368, 618,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2163, 3, 378, 633,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2208, 3, 388, 648,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2253, 3, 398, 663,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2298, 3, 408, 678,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2343, 3, 418, 693,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2388, 3, 7, 8,
                                                                       708, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2394, 3, 8, 9,
                                                                       711, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2400, 3, 9, 10,
                                                                       714, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2406, 3, 10, 11,
                                                                       717, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2412, 3, 11, 12,
                                                                       720, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2418, 3, 12, 13,
                                                                       723, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2424, 3, 13, 14,
                                                                       726, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2430, 3, 14, 15,
                                                                       729, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2436, 3, 15, 16,
                                                                       732, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2442, 3, 16, 17,
                                                                       735, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2448, 3, 17, 18,
                                                                       738, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2454, 3, 21, 22,
                                                                       741, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2460, 3, 22, 23,
                                                                       744, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2466, 3, 23, 24,
                                                                       747, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2472, 3, 24, 25,
                                                                       750, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2478, 3, 25, 26,
                                                                       753, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2484, 3, 26, 27,
                                                                       756, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2490, 3, 27, 28,
                                                                       759, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2496, 3, 28, 29,
                                                                       762, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2502, 3, 29, 30,
                                                                       765, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2508, 3, 30, 31,
                                                                       768, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2514, 3, 31, 32,
                                                                       771, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2520, 0, 3, 2388,
                                                                       708, 2394, 774, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2538, 0, 3, 2394,
                                                                       711, 2400, 783, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2556, 0, 3, 2400,
                                                                       714, 2406, 792, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2574, 0, 3, 2406,
                                                                       717, 2412, 801, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2592, 0, 3, 2412,
                                                                       720, 2418, 810, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2610, 0, 3, 2418,
                                                                       723, 2424, 819, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2628, 0, 3, 2424,
                                                                       726, 2430, 828, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2646, 0, 3, 2430,
                                                                       729, 2436, 837, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2664, 0, 3, 2436,
                                                                       732, 2442, 846, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2682, 0, 3, 2442,
                                                                       735, 2448, 855, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2700, 0, 3, 2454,
                                                                       741, 2460, 864, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2718, 0, 3, 2460,
                                                                       744, 2466, 873, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2736, 0, 3, 2466,
                                                                       747, 2472, 882, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2754, 0, 3, 2472,
                                                                       750, 2478, 891, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2772, 0, 3, 2478,
                                                                       753, 2484, 900, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2790, 0, 3, 2484,
                                                                       756, 2490, 909, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2808, 0, 3, 2490,
                                                                       759, 2496, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2826, 0, 3, 2496,
                                                                       762, 2502, 927, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2844, 0, 3, 2502,
                                                                       765, 2508, 936, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2862, 0, 3, 2508,
                                                                       768, 2514, 945, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2880, 0, 3, 2520,
                                                                       774, 2538, 106, 112, 954,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2916, 0, 3, 2538,
                                                                       783, 2556, 112, 118, 972,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2952, 0, 3, 2556,
                                                                       792, 2574, 118, 124, 990,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2988, 0, 3, 2574,
                                                                       801, 2592, 124, 130, 1008,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3024, 0, 3, 2592,
                                                                       810, 2610, 130, 136, 1026,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3060, 0, 3, 2610,
                                                                       819, 2628, 136, 142, 1044,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3096, 0, 3, 2628,
                                                                       828, 2646, 142, 148, 1062,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3132, 0, 3, 2646,
                                                                       837, 2664, 148, 154, 1080,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3168, 0, 3, 2664,
                                                                       846, 2682, 154, 160, 1098,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3204, 0, 3, 2700,
                                                                       864, 2718, 172, 178, 1116,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3240, 0, 3, 2718,
                                                                       873, 2736, 178, 184, 1134,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3276, 0, 3, 2736,
                                                                       882, 2754, 184, 190, 1152,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3312, 0, 3, 2754,
                                                                       891, 2772, 190, 196, 1170,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3348, 0, 3, 2772,
                                                                       900, 2790, 196, 202, 1188,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3384, 0, 3, 2790,
                                                                       909, 2808, 202, 208, 1206,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3420, 0, 3, 2808,
                                                                       918, 2826, 208, 214, 1224,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3456, 0, 3, 2826,
                                                                       927, 2844, 214, 220, 1242,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3492, 0, 3, 2844,
                                                                       936, 2862, 220, 226, 1260,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3528, 0, 3, 2880,
                                                                       954, 2916, 238, 248, 1278,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3588, 0, 3, 2916,
                                                                       972, 2952, 248, 258, 1308,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3648, 0, 3, 2952,
                                                                       990, 2988, 258, 268, 1338,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3708, 0, 3, 2988,
                                                                       1008, 3024, 268, 278,
                                                                       1368, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3768, 0, 3, 3024,
                                                                       1026, 3060, 278, 288,
                                                                       1398, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3828, 0, 3, 3060,
                                                                       1044, 3096, 288, 298,
                                                                       1428, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3888, 0, 3, 3096,
                                                                       1062, 3132, 298, 308,
                                                                       1458, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3948, 0, 3, 3132,
                                                                       1080, 3168, 308, 318,
                                                                       1488, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4008, 0, 3, 3204,
                                                                       1116, 3240, 338, 348,
                                                                       1518, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4068, 0, 3, 3240,
                                                                       1134, 3276, 348, 358,
                                                                       1548, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4128, 0, 3, 3276,
                                                                       1152, 3312, 358, 368,
                                                                       1578, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4188, 0, 3, 3312,
                                                                       1170, 3348, 368, 378,
                                                                       1608, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4248, 0, 3, 3348,
                                                                       1188, 3384, 378, 388,
                                                                       1638, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4308, 0, 3, 3384,
                                                                       1206, 3420, 388, 398,
                                                                       1668, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4368, 0, 3, 3420,
                                                                       1224, 3456, 398, 408,
                                                                       1698, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4428, 0, 3, 3456,
                                                                       1242, 3492, 408, 418,
                                                                       1728, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4488, 0, 3, 3528,
                                                                       1278, 3588, 438, 453,
                                                                       1758, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4578, 0, 3, 3588,
                                                                       1308, 3648, 453, 468,
                                                                       1803, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4668, 0, 3, 3648,
                                                                       1338, 3708, 468, 483,
                                                                       1848, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4758, 0, 3, 3708,
                                                                       1368, 3768, 483, 498,
                                                                       1893, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4848, 0, 3, 3768,
                                                                       1398, 3828, 498, 513,
                                                                       1938, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4938, 0, 3, 3828,
                                                                       1428, 3888, 513, 528,
                                                                       1983, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5028, 0, 3, 3888,
                                                                       1458, 3948, 528, 543,
                                                                       2028, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5118, 0, 3, 4008,
                                                                       1518, 4068, 573, 588,
                                                                       2073, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5208, 0, 3, 4068,
                                                                       1548, 4128, 588, 603,
                                                                       2118, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5298, 0, 3, 4128,
                                                                       1578, 4188, 603, 618,
                                                                       2163, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5388, 0, 3, 4188,
                                                                       1608, 4248, 618, 633,
                                                                       2208, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5478, 0, 3, 4248,
                                                                       1638, 4308, 633, 648,
                                                                       2253, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5568, 0, 3, 4308,
                                                                       1668, 4368, 648, 663,
                                                                       2298, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5658, 0, 3, 4368,
                                                                       1698, 4428, 663, 678,
                                                                       2343, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5748, 3, 708, 711,
                                                                       2400, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5758, 3, 711, 714,
                                                                       2406, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5768, 3, 714, 717,
                                                                       2412, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5778, 3, 717, 720,
                                                                       2418, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5788, 3, 720, 723,
                                                                       2424, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5798, 3, 723, 726,
                                                                       2430, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5808, 3, 726, 729,
                                                                       2436, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5818, 3, 729, 732,
                                                                       2442, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5828, 3, 732, 735,
                                                                       2448, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5838, 3, 741, 744,
                                                                       2466, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5848, 3, 744, 747,
                                                                       2472, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5858, 3, 747, 750,
                                                                       2478, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5868, 3, 750, 753,
                                                                       2484, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5878, 3, 753, 756,
                                                                       2490, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5888, 3, 756, 759,
                                                                       2496, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5898, 3, 759, 762,
                                                                       2502, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5908, 3, 762, 765,
                                                                       2508, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5918, 3, 765, 768,
                                                                       2514, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5928, 0, 3, 5748,
                                                                       2400, 5758, 774, 783,
                                                                       2556, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5958, 0, 3, 5758,
                                                                       2406, 5768, 783, 792,
                                                                       2574, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5988, 0, 3, 5768,
                                                                       2412, 5778, 792, 801,
                                                                       2592, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6018, 0, 3, 5778,
                                                                       2418, 5788, 801, 810,
                                                                       2610, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6048, 0, 3, 5788,
                                                                       2424, 5798, 810, 819,
                                                                       2628, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6078, 0, 3, 5798,
                                                                       2430, 5808, 819, 828,
                                                                       2646, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6108, 0, 3, 5808,
                                                                       2436, 5818, 828, 837,
                                                                       2664, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6138, 0, 3, 5818,
                                                                       2442, 5828, 837, 846,
                                                                       2682, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6168, 0, 3, 5838,
                                                                       2466, 5848, 864, 873,
                                                                       2736, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6198, 0, 3, 5848,
                                                                       2472, 5858, 873, 882,
                                                                       2754, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6228, 0, 3, 5858,
                                                                       2478, 5868, 882, 891,
                                                                       2772, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6258, 0, 3, 5868,
                                                                       2484, 5878, 891, 900,
                                                                       2790, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6288, 0, 3, 5878,
                                                                       2490, 5888, 900, 909,
                                                                       2808, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6318, 0, 3, 5888,
                                                                       2496, 5898, 909, 918,
                                                                       2826, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6348, 0, 3, 5898,
                                                                       2502, 5908, 918, 927,
                                                                       2844, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6378, 0, 3, 5908,
                                                                       2508, 5918, 927, 936,
                                                                       2862, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6408, 0, 3, 5928,
                                                                       2556, 5958, 954, 972,
                                                                       2952, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6468, 0, 3, 5958,
                                                                       2574, 5988, 972, 990,
                                                                       2988, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6528, 0, 3, 5988,
                                                                       2592, 6018, 990, 1008,
                                                                       3024, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6588, 0, 3, 6018,
                                                                       2610, 6048, 1008, 1026,
                                                                       3060, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6648, 0, 3, 6048,
                                                                       2628, 6078, 1026, 1044,
                                                                       3096, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6708, 0, 3, 6078,
                                                                       2646, 6108, 1044, 1062,
                                                                       3132, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6768, 0, 3, 6108,
                                                                       2664, 6138, 1062, 1080,
                                                                       3168, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6828, 0, 3, 6168,
                                                                       2736, 6198, 1116, 1134,
                                                                       3276, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6888, 0, 3, 6198,
                                                                       2754, 6228, 1134, 1152,
                                                                       3312, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6948, 0, 3, 6228,
                                                                       2772, 6258, 1152, 1170,
                                                                       3348, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7008, 0, 3, 6258,
                                                                       2790, 6288, 1170, 1188,
                                                                       3384, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7068, 0, 3, 6288,
                                                                       2808, 6318, 1188, 1206,
                                                                       3420, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7128, 0, 3, 6318,
                                                                       2826, 6348, 1206, 1224,
                                                                       3456, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7188, 0, 3, 6348,
                                                                       2844, 6378, 1224, 1242,
                                                                       3492, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7248, 0, 3, 6408,
                                                                       2952, 6468, 1278, 1308,
                                                                       3648, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7348, 0, 3, 6468,
                                                                       2988, 6528, 1308, 1338,
                                                                       3708, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7448, 0, 3, 6528,
                                                                       3024, 6588, 1338, 1368,
                                                                       3768, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7548, 0, 3, 6588,
                                                                       3060, 6648, 1368, 1398,
                                                                       3828, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7648, 0, 3, 6648,
                                                                       3096, 6708, 1398, 1428,
                                                                       3888, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7748, 0, 3, 6708,
                                                                       3132, 6768, 1428, 1458,
                                                                       3948, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7848, 0, 3, 6828,
                                                                       3276, 6888, 1518, 1548,
                                                                       4128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7948, 0, 3, 6888,
                                                                       3312, 6948, 1548, 1578,
                                                                       4188, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8048, 0, 3, 6948,
                                                                       3348, 7008, 1578, 1608,
                                                                       4248, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8148, 0, 3, 7008,
                                                                       3384, 7068, 1608, 1638,
                                                                       4308, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8248, 0, 3, 7068,
                                                                       3420, 7128, 1638, 1668,
                                                                       4368, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8348, 0, 3, 7128,
                                                                       3456, 7188, 1668, 1698,
                                                                       4428, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8448, 0, 3, 7248,
                                                                       3648, 7348, 1758, 1803,
                                                                       4668, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8598, 0, 3, 7348,
                                                                       3708, 7448, 1803, 1848,
                                                                       4758, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8748, 0, 3, 7448,
                                                                       3768, 7548, 1848, 1893,
                                                                       4848, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8898, 0, 3, 7548,
                                                                       3828, 7648, 1893, 1938,
                                                                       4938, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9048, 0, 3, 7648,
                                                                       3888, 7748, 1938, 1983,
                                                                       5028, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9198, 0, 3, 7848,
                                                                       4128, 7948, 2073, 2118,
                                                                       5298, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9348, 0, 3, 7948,
                                                                       4188, 8048, 2118, 2163,
                                                                       5388, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9498, 0, 3, 8048,
                                                                       4248, 8148, 2163, 2208,
                                                                       5478, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9648, 0, 3, 8148,
                                                                       4308, 8248, 2208, 2253,
                                                                       5568, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 9798, 0, 3, 8248,
                                                                       4368, 8348, 2253, 2298,
                                                                       5658, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9948, 3, 2388,
                                                                       2394, 5748, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9963, 3, 2394,
                                                                       2400, 5758, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9978, 3, 2400,
                                                                       2406, 5768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9993, 3, 2406,
                                                                       2412, 5778, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10008, 3, 2412,
                                                                       2418, 5788, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10023, 3, 2418,
                                                                       2424, 5798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10038, 3, 2424,
                                                                       2430, 5808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10053, 3, 2430,
                                                                       2436, 5818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10068, 3, 2436,
                                                                       2442, 5828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10083, 3, 2454,
                                                                       2460, 5838, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10098, 3, 2460,
                                                                       2466, 5848, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10113, 3, 2466,
                                                                       2472, 5858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10128, 3, 2472,
                                                                       2478, 5868, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10143, 3, 2478,
                                                                       2484, 5878, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10158, 3, 2484,
                                                                       2490, 5888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10173, 3, 2490,
                                                                       2496, 5898, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10188, 3, 2496,
                                                                       2502, 5908, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10203, 3, 2502,
                                                                       2508, 5918, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10218, 0, 3, 9948,
                                                                       5748, 9963, 2520, 2538,
                                                                       5928, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10263, 0, 3, 9963,
                                                                       5758, 9978, 2538, 2556,
                                                                       5958, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10308, 0, 3, 9978,
                                                                       5768, 9993, 2556, 2574,
                                                                       5988, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10353, 0, 3, 9993,
                                                                       5778, 10008, 2574, 2592,
                                                                       6018, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10398, 0, 3,
                                                                       10008, 5788, 10023, 2592,
                                                                       2610, 6048, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10443, 0, 3,
                                                                       10023, 5798, 10038, 2610,
                                                                       2628, 6078, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10488, 0, 3,
                                                                       10038, 5808, 10053, 2628,
                                                                       2646, 6108, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10533, 0, 3,
                                                                       10053, 5818, 10068, 2646,
                                                                       2664, 6138, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10578, 0, 3,
                                                                       10083, 5838, 10098, 2700,
                                                                       2718, 6168, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10623, 0, 3,
                                                                       10098, 5848, 10113, 2718,
                                                                       2736, 6198, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10668, 0, 3,
                                                                       10113, 5858, 10128, 2736,
                                                                       2754, 6228, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10713, 0, 3,
                                                                       10128, 5868, 10143, 2754,
                                                                       2772, 6258, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10758, 0, 3,
                                                                       10143, 5878, 10158, 2772,
                                                                       2790, 6288, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10803, 0, 3,
                                                                       10158, 5888, 10173, 2790,
                                                                       2808, 6318, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10848, 0, 3,
                                                                       10173, 5898, 10188, 2808,
                                                                       2826, 6348, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10893, 0, 3,
                                                                       10188, 5908, 10203, 2826,
                                                                       2844, 6378, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10938, 0, 3,
                                                                       10218, 5928, 10263, 2880,
                                                                       2916, 6408, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11028, 0, 3,
                                                                       10263, 5958, 10308, 2916,
                                                                       2952, 6468, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11118, 0, 3,
                                                                       10308, 5988, 10353, 2952,
                                                                       2988, 6528, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11208, 0, 3,
                                                                       10353, 6018, 10398, 2988,
                                                                       3024, 6588, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11298, 0, 3,
                                                                       10398, 6048, 10443, 3024,
                                                                       3060, 6648, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11388, 0, 3,
                                                                       10443, 6078, 10488, 3060,
                                                                       3096, 6708, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11478, 0, 3,
                                                                       10488, 6108, 10533, 3096,
                                                                       3132, 6768, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11568, 0, 3,
                                                                       10578, 6168, 10623, 3204,
                                                                       3240, 6828, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11658, 0, 3,
                                                                       10623, 6198, 10668, 3240,
                                                                       3276, 6888, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11748, 0, 3,
                                                                       10668, 6228, 10713, 3276,
                                                                       3312, 6948, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11838, 0, 3,
                                                                       10713, 6258, 10758, 3312,
                                                                       3348, 7008, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11928, 0, 3,
                                                                       10758, 6288, 10803, 3348,
                                                                       3384, 7068, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12018, 0, 3,
                                                                       10803, 6318, 10848, 3384,
                                                                       3420, 7128, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12108, 0, 3,
                                                                       10848, 6348, 10893, 3420,
                                                                       3456, 7188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 12198, 0, 3,
                                                                       10938, 6408, 11028, 3528,
                                                                       3588, 7248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 12348, 0, 3,
                                                                       11028, 6468, 11118, 3588,
                                                                       3648, 7348, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 12498, 0, 3,
                                                                       11118, 6528, 11208, 3648,
                                                                       3708, 7448, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 12648, 0, 3,
                                                                       11208, 6588, 11298, 3708,
                                                                       3768, 7548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 12798, 0, 3,
                                                                       11298, 6648, 11388, 3768,
                                                                       3828, 7648, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 12948, 0, 3,
                                                                       11388, 6708, 11478, 3828,
                                                                       3888, 7748, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13098, 0, 3,
                                                                       11568, 6828, 11658, 4008,
                                                                       4068, 7848, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13248, 0, 3,
                                                                       11658, 6888, 11748, 4068,
                                                                       4128, 7948, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13398, 0, 3,
                                                                       11748, 6948, 11838, 4128,
                                                                       4188, 8048, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13548, 0, 3,
                                                                       11838, 7008, 11928, 4188,
                                                                       4248, 8148, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13698, 0, 3,
                                                                       11928, 7068, 12018, 4248,
                                                                       4308, 8248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13848, 0, 3,
                                                                       12018, 7128, 12108, 4308,
                                                                       4368, 8348, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 13998, 0, 3,
                                                                       12198, 7248, 12348, 4488,
                                                                       4578, 8448, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14223, 0, 3,
                                                                       12348, 7348, 12498, 4578,
                                                                       4668, 8598, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14448, 0, 3,
                                                                       12498, 7448, 12648, 4668,
                                                                       4758, 8748, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14673, 0, 3,
                                                                       12648, 7548, 12798, 4758,
                                                                       4848, 8898, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 14898, 0, 3,
                                                                       12798, 7648, 12948, 4848,
                                                                       4938, 9048, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 15123, 0, 3,
                                                                       13098, 7848, 13248, 5118,
                                                                       5208, 9198, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 15348, 0, 3,
                                                                       13248, 7948, 13398, 5208,
                                                                       5298, 9348, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 15573, 0, 3,
                                                                       13398, 8048, 13548, 5298,
                                                                       5388, 9498, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 15798, 0, 3,
                                                                       13548, 8148, 13698, 5388,
                                                                       5478, 9648, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 16023, 0, 3,
                                                                       13698, 8248, 13848, 5478,
                                                                       5568, 9798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16248, 3, 5748,
                                                                       5758, 9978, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16269, 3, 5758,
                                                                       5768, 9993, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16290, 3, 5768,
                                                                       5778, 10008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16311, 3, 5778,
                                                                       5788, 10023, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16332, 3, 5788,
                                                                       5798, 10038, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16353, 3, 5798,
                                                                       5808, 10053, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16374, 3, 5808,
                                                                       5818, 10068, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16395, 3, 5838,
                                                                       5848, 10113, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16416, 3, 5848,
                                                                       5858, 10128, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16437, 3, 5858,
                                                                       5868, 10143, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16458, 3, 5868,
                                                                       5878, 10158, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16479, 3, 5878,
                                                                       5888, 10173, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16500, 3, 5888,
                                                                       5898, 10188, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 16521, 3, 5898,
                                                                       5908, 10203, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 16542, 0, 3,
                                                                       16248, 9978, 16269, 5928,
                                                                       5958, 10308, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 16605, 0, 3,
                                                                       16269, 9993, 16290, 5958,
                                                                       5988, 10353, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 16668, 0, 3,
                                                                       16290, 10008, 16311, 5988,
                                                                       6018, 10398, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 16731, 0, 3,
                                                                       16311, 10023, 16332, 6018,
                                                                       6048, 10443, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 16794, 0, 3,
                                                                       16332, 10038, 16353, 6048,
                                                                       6078, 10488, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 16857, 0, 3,
                                                                       16353, 10053, 16374, 6078,
                                                                       6108, 10533, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 16920, 0, 3,
                                                                       16395, 10113, 16416, 6168,
                                                                       6198, 10668, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 16983, 0, 3,
                                                                       16416, 10128, 16437, 6198,
                                                                       6228, 10713, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17046, 0, 3,
                                                                       16437, 10143, 16458, 6228,
                                                                       6258, 10758, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17109, 0, 3,
                                                                       16458, 10158, 16479, 6258,
                                                                       6288, 10803, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17172, 0, 3,
                                                                       16479, 10173, 16500, 6288,
                                                                       6318, 10848, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 17235, 0, 3,
                                                                       16500, 10188, 16521, 6318,
                                                                       6348, 10893, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 17298, 0, 3,
                                                                       16542, 10308, 16605, 6408,
                                                                       6468, 11118, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 17424, 0, 3,
                                                                       16605, 10353, 16668, 6468,
                                                                       6528, 11208, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 17550, 0, 3,
                                                                       16668, 10398, 16731, 6528,
                                                                       6588, 11298, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 17676, 0, 3,
                                                                       16731, 10443, 16794, 6588,
                                                                       6648, 11388, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 17802, 0, 3,
                                                                       16794, 10488, 16857, 6648,
                                                                       6708, 11478, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 17928, 0, 3,
                                                                       16920, 10668, 16983, 6828,
                                                                       6888, 11748, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18054, 0, 3,
                                                                       16983, 10713, 17046, 6888,
                                                                       6948, 11838, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18180, 0, 3,
                                                                       17046, 10758, 17109, 6948,
                                                                       7008, 11928, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18306, 0, 3,
                                                                       17109, 10803, 17172, 7008,
                                                                       7068, 12018, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 18432, 0, 3,
                                                                       17172, 10848, 17235, 7068,
                                                                       7128, 12108, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 18558, 0, 3,
                                                                       17298, 11118, 17424, 7248,
                                                                       7348, 12498, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 18768, 0, 3,
                                                                       17424, 11208, 17550, 7348,
                                                                       7448, 12648, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 18978, 0, 3,
                                                                       17550, 11298, 17676, 7448,
                                                                       7548, 12798, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19188, 0, 3,
                                                                       17676, 11388, 17802, 7548,
                                                                       7648, 12948, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19398, 0, 3,
                                                                       17928, 11748, 18054, 7848,
                                                                       7948, 13398, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19608, 0, 3,
                                                                       18054, 11838, 18180, 7948,
                                                                       8048, 13548, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 19818, 0, 3,
                                                                       18180, 11928, 18306, 8048,
                                                                       8148, 13698, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 20028, 0, 3,
                                                                       18306, 12018, 18432, 8148,
                                                                       8248, 13848, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 20238, 0, 3,
                                                                       18558, 12498, 18768, 8448,
                                                                       8598, 14448, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 20553, 0, 3,
                                                                       18768, 12648, 18978, 8598,
                                                                       8748, 14673, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 20868, 0, 3,
                                                                       18978, 12798, 19188, 8748,
                                                                       8898, 14898, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 21183, 0, 3,
                                                                       19398, 13398, 19608, 9198,
                                                                       9348, 15573, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 21498, 0, 3,
                                                                       19608, 13548, 19818, 9348,
                                                                       9498, 15798, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 21813, 0, 3,
                                                                       19818, 13698, 20028, 9498,
                                                                       9648, 16023, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22128, 3, 9948,
                                                                       9963, 16248, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22156, 3, 9963,
                                                                       9978, 16269, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22184, 3, 9978,
                                                                       9993, 16290, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22212, 3, 9993,
                                                                       10008, 16311, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22240, 3, 10008,
                                                                       10023, 16332, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22268, 3, 10023,
                                                                       10038, 16353, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22296, 3, 10038,
                                                                       10053, 16374, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22324, 3, 10083,
                                                                       10098, 16395, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22352, 3, 10098,
                                                                       10113, 16416, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22380, 3, 10113,
                                                                       10128, 16437, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22408, 3, 10128,
                                                                       10143, 16458, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22436, 3, 10143,
                                                                       10158, 16479, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22464, 3, 10158,
                                                                       10173, 16500, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 22492, 3, 10173,
                                                                       10188, 16521, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 22520, 0, 3,
                                                                       22128, 16248, 22156,
                                                                       10218, 10263, 16542,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 22604, 0, 3,
                                                                       22156, 16269, 22184,
                                                                       10263, 10308, 16605,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 22688, 0, 3,
                                                                       22184, 16290, 22212,
                                                                       10308, 10353, 16668,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 22772, 0, 3,
                                                                       22212, 16311, 22240,
                                                                       10353, 10398, 16731,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 22856, 0, 3,
                                                                       22240, 16332, 22268,
                                                                       10398, 10443, 16794,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 22940, 0, 3,
                                                                       22268, 16353, 22296,
                                                                       10443, 10488, 16857,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 23024, 0, 3,
                                                                       22324, 16395, 22352,
                                                                       10578, 10623, 16920,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 23108, 0, 3,
                                                                       22352, 16416, 22380,
                                                                       10623, 10668, 16983,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 23192, 0, 3,
                                                                       22380, 16437, 22408,
                                                                       10668, 10713, 17046,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 23276, 0, 3,
                                                                       22408, 16458, 22436,
                                                                       10713, 10758, 17109,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 23360, 0, 3,
                                                                       22436, 16479, 22464,
                                                                       10758, 10803, 17172,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 23444, 0, 3,
                                                                       22464, 16500, 22492,
                                                                       10803, 10848, 17235,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 23528, 0, 3,
                                                                       22520, 16542, 22604,
                                                                       10938, 11028, 17298,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 23696, 0, 3,
                                                                       22604, 16605, 22688,
                                                                       11028, 11118, 17424,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 23864, 0, 3,
                                                                       22688, 16668, 22772,
                                                                       11118, 11208, 17550,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 24032, 0, 3,
                                                                       22772, 16731, 22856,
                                                                       11208, 11298, 17676,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 24200, 0, 3,
                                                                       22856, 16794, 22940,
                                                                       11298, 11388, 17802,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 24368, 0, 3,
                                                                       23024, 16920, 23108,
                                                                       11568, 11658, 17928,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 24536, 0, 3,
                                                                       23108, 16983, 23192,
                                                                       11658, 11748, 18054,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 24704, 0, 3,
                                                                       23192, 17046, 23276,
                                                                       11748, 11838, 18180,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 24872, 0, 3,
                                                                       23276, 17109, 23360,
                                                                       11838, 11928, 18306,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 25040, 0, 3,
                                                                       23360, 17172, 23444,
                                                                       11928, 12018, 18432,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 25208, 0, 3,
                                                                       23528, 17298, 23696,
                                                                       12198, 12348, 18558,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 25488, 0, 3,
                                                                       23696, 17424, 23864,
                                                                       12348, 12498, 18768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 25768, 0, 3,
                                                                       23864, 17550, 24032,
                                                                       12498, 12648, 18978,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 26048, 0, 3,
                                                                       24032, 17676, 24200,
                                                                       12648, 12798, 19188,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 26328, 0, 3,
                                                                       24368, 17928, 24536,
                                                                       13098, 13248, 19398,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 26608, 0, 3,
                                                                       24536, 18054, 24704,
                                                                       13248, 13398, 19608,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 26888, 0, 3,
                                                                       24704, 18180, 24872,
                                                                       13398, 13548, 19818,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 27168, 0, 3,
                                                                       24872, 18306, 25040,
                                                                       13548, 13698, 20028,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 27448, 0, 3,
                                                                       25208, 18558, 25488,
                                                                       13998, 14223, 20238,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 27868, 0, 3,
                                                                       25488, 18768, 25768,
                                                                       14223, 14448, 20553,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 28288, 0, 3,
                                                                       25768, 18978, 26048,
                                                                       14448, 14673, 20868,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 28708, 0, 3,
                                                                       26328, 19398, 26608,
                                                                       15123, 15348, 21183,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 29128, 0, 3,
                                                                       26608, 19608, 26888,
                                                                       15348, 15573, 21498,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 29548, 0, 3,
                                                                       26888, 19818, 27168,
                                                                       15573, 15798, 21813,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 29968, 3, 16248,
                                                                       16269, 22184, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 30004, 3, 16269,
                                                                       16290, 22212, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 30040, 3, 16290,
                                                                       16311, 22240, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 30076, 3, 16311,
                                                                       16332, 22268, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 30112, 3, 16332,
                                                                       16353, 22296, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 30148, 3, 16395,
                                                                       16416, 22380, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 30184, 3, 16416,
                                                                       16437, 22408, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 30220, 3, 16437,
                                                                       16458, 22436, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 30256, 3, 16458,
                                                                       16479, 22464, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 30292, 3, 16479,
                                                                       16500, 22492, ncols,
                                                                       gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 30328, 0, 3,
                                                                       29968, 22184, 30004,
                                                                       16542, 16605, 22688,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 30436, 0, 3,
                                                                       30004, 22212, 30040,
                                                                       16605, 16668, 22772,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 30544, 0, 3,
                                                                       30040, 22240, 30076,
                                                                       16668, 16731, 22856,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 30652, 0, 3,
                                                                       30076, 22268, 30112,
                                                                       16731, 16794, 22940,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 30760, 0, 3,
                                                                       30148, 22380, 30184,
                                                                       16920, 16983, 23192,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 30868, 0, 3,
                                                                       30184, 22408, 30220,
                                                                       16983, 17046, 23276,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 30976, 0, 3,
                                                                       30220, 22436, 30256,
                                                                       17046, 17109, 23360,
                                                                       ncols, gamma, p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 31084, 0, 3,
                                                                       30256, 22464, 30292,
                                                                       17109, 17172, 23444,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 31192, 0, 3,
                                                                       30328, 22688, 30436,
                                                                       17298, 17424, 23864,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 31408, 0, 3,
                                                                       30436, 22772, 30544,
                                                                       17424, 17550, 24032,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 31624, 0, 3,
                                                                       30544, 22856, 30652,
                                                                       17550, 17676, 24200,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 31840, 0, 3,
                                                                       30760, 23192, 30868,
                                                                       17928, 18054, 24704,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 32056, 0, 3,
                                                                       30868, 23276, 30976,
                                                                       18054, 18180, 24872,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 32272, 0, 3,
                                                                       30976, 23360, 31084,
                                                                       18180, 18306, 25040,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 32488, 0, 3,
                                                                       31192, 23864, 31408,
                                                                       18558, 18768, 25768,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 32848, 0, 3,
                                                                       31408, 24032, 31624,
                                                                       18768, 18978, 26048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 33208, 0, 3,
                                                                       31840, 24704, 32056,
                                                                       19398, 19608, 26888,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 33568, 0, 3,
                                                                       32056, 24872, 32272,
                                                                       19608, 19818, 27168,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 33928, 0, 3,
                                                                       32488, 25768, 32848,
                                                                       20238, 20553, 28288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 34468, 0, 3,
                                                                       33208, 26888, 33568,
                                                                       21183, 21498, 29548,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 35008, 3, 22128,
                                                                       22156, 29968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 35053, 3, 22156,
                                                                       22184, 30004, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 35098, 3, 22184,
                                                                       22212, 30040, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 35143, 3, 22212,
                                                                       22240, 30076, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 35188, 3, 22240,
                                                                       22268, 30112, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 35233, 3, 22324,
                                                                       22352, 30148, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 35278, 3, 22352,
                                                                       22380, 30184, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 35323, 3, 22380,
                                                                       22408, 30220, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 35368, 3, 22408,
                                                                       22436, 30256, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 35413, 3, 22436,
                                                                       22464, 30292, ncols,
                                                                       gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 35458, 0, 3,
                                                                       35008, 29968, 35053,
                                                                       22520, 22604, 30328,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 35593, 0, 3,
                                                                       35053, 30004, 35098,
                                                                       22604, 22688, 30436,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 35728, 0, 3,
                                                                       35098, 30040, 35143,
                                                                       22688, 22772, 30544,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 35863, 0, 3,
                                                                       35143, 30076, 35188,
                                                                       22772, 22856, 30652,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 35998, 0, 3,
                                                                       35233, 30148, 35278,
                                                                       23024, 23108, 30760,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 36133, 0, 3,
                                                                       35278, 30184, 35323,
                                                                       23108, 23192, 30868,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 36268, 0, 3,
                                                                       35323, 30220, 35368,
                                                                       23192, 23276, 30976,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 36403, 0, 3,
                                                                       35368, 30256, 35413,
                                                                       23276, 23360, 31084,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 36538, 0, 3,
                                                                       35458, 30328, 35593,
                                                                       23528, 23696, 31192,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 36808, 0, 3,
                                                                       35593, 30436, 35728,
                                                                       23696, 23864, 31408,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 37078, 0, 3,
                                                                       35728, 30544, 35863,
                                                                       23864, 24032, 31624,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 37348, 0, 3,
                                                                       35998, 30760, 36133,
                                                                       24368, 24536, 31840,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 37618, 0, 3,
                                                                       36133, 30868, 36268,
                                                                       24536, 24704, 32056,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 37888, 0, 3,
                                                                       36268, 30976, 36403,
                                                                       24704, 24872, 32272,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 38158, 0, 3,
                                                                       36538, 31192, 36808,
                                                                       25208, 25488, 32488,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 38608, 0, 3,
                                                                       36808, 31408, 37078,
                                                                       25488, 25768, 32848,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 39058, 0, 3,
                                                                       37348, 31840, 37618,
                                                                       26328, 26608, 33208,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 39508, 0, 3,
                                                                       37618, 32056, 37888,
                                                                       26608, 26888, 33568,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 39958, 0, 3,
                                                                       38158, 32488, 38608,
                                                                       27448, 27868, 33928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 40633, 0, 3,
                                                                       39058, 33208, 39508,
                                                                       28708, 29128, 34468,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 41308, 36538, 270, ncols);

                    simdfunc::contract_primitives(buffer, 41680, 37348, 270, ncols);

                    simdfunc::contract_primitives(buffer, 42052, 38158, 450, ncols);

                    simdfunc::contract_primitives(buffer, 42672, 39058, 450, ncols);

                    simdfunc::contract_primitives(buffer, 43292, 39958, 675, ncols);

                    simdfunc::contract_primitives(buffer, 44222, 40633, 675, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 41578, 41308, 6, 1, nmax);

        simdtrf::transform_l_inner(buffer, 41950, 41680, 6, 1, nmax);

        simdtrf::transform_l_inner(buffer, 42502, 42052, 10, 1, nmax);

        simdtrf::transform_l_inner(buffer, 43122, 42672, 10, 1, nmax);

        simdtrf::transform_l_inner(buffer, 43967, 43292, 15, 1, nmax);

        simdtrf::transform_l_inner(buffer, 44897, 44222, 15, 1, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 45152, 41578, 42502, 17,
                                             nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 45458, 41950, 43122, 17,
                                             nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 45764, 42502, 43967, 17,
                                             nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 46274, 43122, 44897, 17,
                                             nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 46784, 45152, 45764, 17,
                                             nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 47396, 45458, 46274, 17,
                                             nmax);

        simdtrf::transform_d_inner(buffer, 48008, 47396, 6, 17, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 48008, 85, nmax);

        simdtrf::transform_d_inner(buffer, 48008, 46784, 6, 17, nmax);

        simdtrf::transform_d_outer(values + 425 * nvalues + n * npairs, nvalues, buffer, 48008,
                                   85, nmax);
    }

    for (size_t m = 0; m < 850; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
