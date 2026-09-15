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


#include "SimdThreeCenterElectronRepulsionRecDDL.hpp"

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
compute_ddl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_ddl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 24518, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 425 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 24518, 20658, 1667, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 12,
                                                             ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 21, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 24, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 27, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 8, 9,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 63, 0, 3, 9, 10,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 69, 0, 3, 10, 11,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 75, 0, 3, 11, 12,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 81, 0, 3, 12, 13,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 87, 0, 3, 13, 14,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 93, 0, 3, 14, 15,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 99, 0, 3, 15, 16,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 105, 0, 3, 16, 17,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 111, 0, 3, 17, 18,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 117, 0, 3, 18, 19,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 123, 0, 3, 21, 24,
                                                                       57, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 133, 0, 3, 24, 27,
                                                                       63, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 27, 30,
                                                                       69, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 153, 0, 3, 30, 33,
                                                                       75, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 163, 0, 3, 33, 36,
                                                                       81, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 36, 39,
                                                                       87, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 183, 0, 3, 39, 42,
                                                                       93, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 42, 45,
                                                                       99, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 45, 48,
                                                                       105, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 48, 51,
                                                                       111, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 223, 0, 3, 57, 63,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 63, 69,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 253, 0, 3, 69, 75,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 75, 81,
                                                                       153, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 283, 0, 3, 81, 87,
                                                                       163, 173, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 87, 93,
                                                                       173, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 313, 0, 3, 93, 99,
                                                                       183, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 99,
                                                                       105, 193, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 343, 0, 3, 105,
                                                                       111, 203, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 358, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 361, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 364, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 367, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 370, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 373, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 376, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 379, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 382, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 385, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 388, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 391, 3, 10, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 400, 3, 11, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 409, 3, 12, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 418, 3, 13, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 427, 3, 14, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 436, 3, 15, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 445, 3, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 454, 3, 17, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 463, 3, 18, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 472, 3, 19, 54,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 481, 3, 27, 69,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 499, 3, 30, 75,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 517, 3, 33, 81,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 535, 3, 36, 87,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 553, 3, 39, 93,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 571, 3, 42, 99,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 589, 3, 45, 105,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 607, 3, 48, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 625, 3, 51, 117,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 643, 3, 69, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 673, 3, 75, 153,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 703, 3, 81, 163,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 733, 3, 87, 173,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 763, 3, 93, 183,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 793, 3, 99, 193,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 823, 3, 105, 203,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 853, 3, 111, 213,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 883, 3, 143, 253,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 928, 3, 153, 268,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 973, 3, 163, 283,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1018, 3, 173, 298,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1063, 3, 183, 313,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1108, 3, 193, 328,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1153, 3, 203, 343,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1198, 3, 8, 9,
                                                                       358, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1204, 3, 9, 10,
                                                                       361, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1210, 3, 10, 11,
                                                                       364, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1216, 3, 11, 12,
                                                                       367, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1222, 3, 12, 13,
                                                                       370, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1228, 3, 13, 14,
                                                                       373, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1234, 3, 14, 15,
                                                                       376, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1240, 3, 15, 16,
                                                                       379, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1246, 3, 16, 17,
                                                                       382, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1252, 3, 17, 18,
                                                                       385, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1258, 3, 18, 19,
                                                                       388, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1264, 0, 3, 1198,
                                                                       358, 1204, 391, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1282, 0, 3, 1204,
                                                                       361, 1210, 400, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1300, 0, 3, 1210,
                                                                       364, 1216, 409, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1318, 0, 3, 1216,
                                                                       367, 1222, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1336, 0, 3, 1222,
                                                                       370, 1228, 427, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1354, 0, 3, 1228,
                                                                       373, 1234, 436, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1372, 0, 3, 1234,
                                                                       376, 1240, 445, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1390, 0, 3, 1240,
                                                                       379, 1246, 454, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 1246,
                                                                       382, 1252, 463, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1426, 0, 3, 1252,
                                                                       385, 1258, 472, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1444, 0, 3, 1264,
                                                                       391, 1282, 57, 63, 481,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1480, 0, 3, 1282,
                                                                       400, 1300, 63, 69, 499,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1516, 0, 3, 1300,
                                                                       409, 1318, 69, 75, 517,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1552, 0, 3, 1318,
                                                                       418, 1336, 75, 81, 535,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1588, 0, 3, 1336,
                                                                       427, 1354, 81, 87, 553,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1624, 0, 3, 1354,
                                                                       436, 1372, 87, 93, 571,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1660, 0, 3, 1372,
                                                                       445, 1390, 93, 99, 589,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1696, 0, 3, 1390,
                                                                       454, 1408, 99, 105, 607,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1732, 0, 3, 1408,
                                                                       463, 1426, 105, 111, 625,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1768, 0, 3, 1444,
                                                                       481, 1480, 123, 133, 643,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1828, 0, 3, 1480,
                                                                       499, 1516, 133, 143, 673,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1888, 0, 3, 1516,
                                                                       517, 1552, 143, 153, 703,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1948, 0, 3, 1552,
                                                                       535, 1588, 153, 163, 733,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2008, 0, 3, 1588,
                                                                       553, 1624, 163, 173, 763,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2068, 0, 3, 1624,
                                                                       571, 1660, 173, 183, 793,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2128, 0, 3, 1660,
                                                                       589, 1696, 183, 193, 823,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2188, 0, 3, 1696,
                                                                       607, 1732, 193, 203, 853,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2248, 0, 3, 1768,
                                                                       643, 1828, 223, 238, 883,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2338, 0, 3, 1828,
                                                                       673, 1888, 238, 253, 928,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2428, 0, 3, 1888,
                                                                       703, 1948, 253, 268, 973,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2518, 0, 3, 1948,
                                                                       733, 2008, 268, 283, 1018,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2608, 0, 3, 2008,
                                                                       763, 2068, 283, 298, 1063,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2698, 0, 3, 2068,
                                                                       793, 2128, 298, 313, 1108,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2788, 0, 3, 2128,
                                                                       823, 2188, 313, 328, 1153,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2878, 3, 358, 361,
                                                                       1210, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2888, 3, 361, 364,
                                                                       1216, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2898, 3, 364, 367,
                                                                       1222, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2908, 3, 367, 370,
                                                                       1228, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2918, 3, 370, 373,
                                                                       1234, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2928, 3, 373, 376,
                                                                       1240, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2938, 3, 376, 379,
                                                                       1246, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2948, 3, 379, 382,
                                                                       1252, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2958, 3, 382, 385,
                                                                       1258, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2968, 0, 3, 2878,
                                                                       1210, 2888, 391, 400,
                                                                       1300, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2998, 0, 3, 2888,
                                                                       1216, 2898, 400, 409,
                                                                       1318, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3028, 0, 3, 2898,
                                                                       1222, 2908, 409, 418,
                                                                       1336, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3058, 0, 3, 2908,
                                                                       1228, 2918, 418, 427,
                                                                       1354, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3088, 0, 3, 2918,
                                                                       1234, 2928, 427, 436,
                                                                       1372, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3118, 0, 3, 2928,
                                                                       1240, 2938, 436, 445,
                                                                       1390, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3148, 0, 3, 2938,
                                                                       1246, 2948, 445, 454,
                                                                       1408, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3178, 0, 3, 2948,
                                                                       1252, 2958, 454, 463,
                                                                       1426, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3208, 0, 3, 2968,
                                                                       1300, 2998, 481, 499,
                                                                       1516, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3268, 0, 3, 2998,
                                                                       1318, 3028, 499, 517,
                                                                       1552, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3328, 0, 3, 3028,
                                                                       1336, 3058, 517, 535,
                                                                       1588, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3388, 0, 3, 3058,
                                                                       1354, 3088, 535, 553,
                                                                       1624, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3448, 0, 3, 3088,
                                                                       1372, 3118, 553, 571,
                                                                       1660, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3508, 0, 3, 3118,
                                                                       1390, 3148, 571, 589,
                                                                       1696, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3568, 0, 3, 3148,
                                                                       1408, 3178, 589, 607,
                                                                       1732, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3628, 0, 3, 3208,
                                                                       1516, 3268, 643, 673,
                                                                       1888, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3728, 0, 3, 3268,
                                                                       1552, 3328, 673, 703,
                                                                       1948, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3828, 0, 3, 3328,
                                                                       1588, 3388, 703, 733,
                                                                       2008, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3928, 0, 3, 3388,
                                                                       1624, 3448, 733, 763,
                                                                       2068, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4028, 0, 3, 3448,
                                                                       1660, 3508, 763, 793,
                                                                       2128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4128, 0, 3, 3508,
                                                                       1696, 3568, 793, 823,
                                                                       2188, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4228, 0, 3, 3628,
                                                                       1888, 3728, 883, 928,
                                                                       2428, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4378, 0, 3, 3728,
                                                                       1948, 3828, 928, 973,
                                                                       2518, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4528, 0, 3, 3828,
                                                                       2008, 3928, 973, 1018,
                                                                       2608, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4678, 0, 3, 3928,
                                                                       2068, 4028, 1018, 1063,
                                                                       2698, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4828, 0, 3, 4028,
                                                                       2128, 4128, 1063, 1108,
                                                                       2788, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4978, 3, 1198,
                                                                       1204, 2878, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4993, 3, 1204,
                                                                       1210, 2888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5008, 3, 1210,
                                                                       1216, 2898, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5023, 3, 1216,
                                                                       1222, 2908, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5038, 3, 1222,
                                                                       1228, 2918, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5053, 3, 1228,
                                                                       1234, 2928, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5068, 3, 1234,
                                                                       1240, 2938, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5083, 3, 1240,
                                                                       1246, 2948, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5098, 3, 1246,
                                                                       1252, 2958, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5113, 0, 3, 4978,
                                                                       2878, 4993, 1264, 1282,
                                                                       2968, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5158, 0, 3, 4993,
                                                                       2888, 5008, 1282, 1300,
                                                                       2998, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5203, 0, 3, 5008,
                                                                       2898, 5023, 1300, 1318,
                                                                       3028, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5248, 0, 3, 5023,
                                                                       2908, 5038, 1318, 1336,
                                                                       3058, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5293, 0, 3, 5038,
                                                                       2918, 5053, 1336, 1354,
                                                                       3088, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5338, 0, 3, 5053,
                                                                       2928, 5068, 1354, 1372,
                                                                       3118, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5383, 0, 3, 5068,
                                                                       2938, 5083, 1372, 1390,
                                                                       3148, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5428, 0, 3, 5083,
                                                                       2948, 5098, 1390, 1408,
                                                                       3178, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5473, 0, 3, 5113,
                                                                       2968, 5158, 1444, 1480,
                                                                       3208, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5563, 0, 3, 5158,
                                                                       2998, 5203, 1480, 1516,
                                                                       3268, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5653, 0, 3, 5203,
                                                                       3028, 5248, 1516, 1552,
                                                                       3328, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5743, 0, 3, 5248,
                                                                       3058, 5293, 1552, 1588,
                                                                       3388, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5833, 0, 3, 5293,
                                                                       3088, 5338, 1588, 1624,
                                                                       3448, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5923, 0, 3, 5338,
                                                                       3118, 5383, 1624, 1660,
                                                                       3508, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 6013, 0, 3, 5383,
                                                                       3148, 5428, 1660, 1696,
                                                                       3568, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6103, 0, 3, 5473,
                                                                       3208, 5563, 1768, 1828,
                                                                       3628, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6253, 0, 3, 5563,
                                                                       3268, 5653, 1828, 1888,
                                                                       3728, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6403, 0, 3, 5653,
                                                                       3328, 5743, 1888, 1948,
                                                                       3828, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6553, 0, 3, 5743,
                                                                       3388, 5833, 1948, 2008,
                                                                       3928, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6703, 0, 3, 5833,
                                                                       3448, 5923, 2008, 2068,
                                                                       4028, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6853, 0, 3, 5923,
                                                                       3508, 6013, 2068, 2128,
                                                                       4128, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 7003, 0, 3, 6103,
                                                                       3628, 6253, 2248, 2338,
                                                                       4228, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 7228, 0, 3, 6253,
                                                                       3728, 6403, 2338, 2428,
                                                                       4378, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 7453, 0, 3, 6403,
                                                                       3828, 6553, 2428, 2518,
                                                                       4528, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 7678, 0, 3, 6553,
                                                                       3928, 6703, 2518, 2608,
                                                                       4678, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 7903, 0, 3, 6703,
                                                                       4028, 6853, 2608, 2698,
                                                                       4828, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8128, 3, 2878,
                                                                       2888, 5008, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8149, 3, 2888,
                                                                       2898, 5023, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8170, 3, 2898,
                                                                       2908, 5038, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8191, 3, 2908,
                                                                       2918, 5053, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8212, 3, 2918,
                                                                       2928, 5068, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8233, 3, 2928,
                                                                       2938, 5083, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8254, 3, 2938,
                                                                       2948, 5098, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 8275, 0, 3, 8128,
                                                                       5008, 8149, 2968, 2998,
                                                                       5203, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 8338, 0, 3, 8149,
                                                                       5023, 8170, 2998, 3028,
                                                                       5248, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 8401, 0, 3, 8170,
                                                                       5038, 8191, 3028, 3058,
                                                                       5293, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 8464, 0, 3, 8191,
                                                                       5053, 8212, 3058, 3088,
                                                                       5338, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 8527, 0, 3, 8212,
                                                                       5068, 8233, 3088, 3118,
                                                                       5383, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 8590, 0, 3, 8233,
                                                                       5083, 8254, 3118, 3148,
                                                                       5428, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 8653, 0, 3, 8275,
                                                                       5203, 8338, 3208, 3268,
                                                                       5653, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 8779, 0, 3, 8338,
                                                                       5248, 8401, 3268, 3328,
                                                                       5743, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 8905, 0, 3, 8401,
                                                                       5293, 8464, 3328, 3388,
                                                                       5833, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 9031, 0, 3, 8464,
                                                                       5338, 8527, 3388, 3448,
                                                                       5923, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 9157, 0, 3, 8527,
                                                                       5383, 8590, 3448, 3508,
                                                                       6013, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 9283, 0, 3, 8653,
                                                                       5653, 8779, 3628, 3728,
                                                                       6403, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 9493, 0, 3, 8779,
                                                                       5743, 8905, 3728, 3828,
                                                                       6553, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 9703, 0, 3, 8905,
                                                                       5833, 9031, 3828, 3928,
                                                                       6703, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 9913, 0, 3, 9031,
                                                                       5923, 9157, 3928, 4028,
                                                                       6853, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 10123, 0, 3, 9283,
                                                                       6403, 9493, 4228, 4378,
                                                                       7453, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 10438, 0, 3, 9493,
                                                                       6553, 9703, 4378, 4528,
                                                                       7678, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 10753, 0, 3, 9703,
                                                                       6703, 9913, 4528, 4678,
                                                                       7903, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 11068, 3, 4978,
                                                                       4993, 8128, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 11096, 3, 4993,
                                                                       5008, 8149, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 11124, 3, 5008,
                                                                       5023, 8170, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 11152, 3, 5023,
                                                                       5038, 8191, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 11180, 3, 5038,
                                                                       5053, 8212, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 11208, 3, 5053,
                                                                       5068, 8233, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 11236, 3, 5068,
                                                                       5083, 8254, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 11264, 0, 3,
                                                                       11068, 8128, 11096, 5113,
                                                                       5158, 8275, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 11348, 0, 3,
                                                                       11096, 8149, 11124, 5158,
                                                                       5203, 8338, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 11432, 0, 3,
                                                                       11124, 8170, 11152, 5203,
                                                                       5248, 8401, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 11516, 0, 3,
                                                                       11152, 8191, 11180, 5248,
                                                                       5293, 8464, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 11600, 0, 3,
                                                                       11180, 8212, 11208, 5293,
                                                                       5338, 8527, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 11684, 0, 3,
                                                                       11208, 8233, 11236, 5338,
                                                                       5383, 8590, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 11768, 0, 3,
                                                                       11264, 8275, 11348, 5473,
                                                                       5563, 8653, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 11936, 0, 3,
                                                                       11348, 8338, 11432, 5563,
                                                                       5653, 8779, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 12104, 0, 3,
                                                                       11432, 8401, 11516, 5653,
                                                                       5743, 8905, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 12272, 0, 3,
                                                                       11516, 8464, 11600, 5743,
                                                                       5833, 9031, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 12440, 0, 3,
                                                                       11600, 8527, 11684, 5833,
                                                                       5923, 9157, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 12608, 0, 3,
                                                                       11768, 8653, 11936, 6103,
                                                                       6253, 9283, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 12888, 0, 3,
                                                                       11936, 8779, 12104, 6253,
                                                                       6403, 9493, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 13168, 0, 3,
                                                                       12104, 8905, 12272, 6403,
                                                                       6553, 9703, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 13448, 0, 3,
                                                                       12272, 9031, 12440, 6553,
                                                                       6703, 9913, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 13728, 0, 3,
                                                                       12608, 9283, 12888, 7003,
                                                                       7228, 10123, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 14148, 0, 3,
                                                                       12888, 9493, 13168, 7228,
                                                                       7453, 10438, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 14568, 0, 3,
                                                                       13168, 9703, 13448, 7453,
                                                                       7678, 10753, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 14988, 3, 8128,
                                                                       8149, 11124, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 15024, 3, 8149,
                                                                       8170, 11152, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 15060, 3, 8170,
                                                                       8191, 11180, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 15096, 3, 8191,
                                                                       8212, 11208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 15132, 3, 8212,
                                                                       8233, 11236, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 15168, 0, 3,
                                                                       14988, 11124, 15024, 8275,
                                                                       8338, 11432, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 15276, 0, 3,
                                                                       15024, 11152, 15060, 8338,
                                                                       8401, 11516, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 15384, 0, 3,
                                                                       15060, 11180, 15096, 8401,
                                                                       8464, 11600, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 15492, 0, 3,
                                                                       15096, 11208, 15132, 8464,
                                                                       8527, 11684, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 15600, 0, 3,
                                                                       15168, 11432, 15276, 8653,
                                                                       8779, 12104, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 15816, 0, 3,
                                                                       15276, 11516, 15384, 8779,
                                                                       8905, 12272, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 16032, 0, 3,
                                                                       15384, 11600, 15492, 8905,
                                                                       9031, 12440, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 16248, 0, 3,
                                                                       15600, 12104, 15816, 9283,
                                                                       9493, 13168, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 16608, 0, 3,
                                                                       15816, 12272, 16032, 9493,
                                                                       9703, 13448, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsk_three_center_electron_repulsion_0(buffer, 16968, 0, 3,
                                                                       16248, 13168, 16608,
                                                                       10123, 10438, 14568,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 17508, 3, 11068,
                                                                       11096, 14988, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 17553, 3, 11096,
                                                                       11124, 15024, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 17598, 3, 11124,
                                                                       11152, 15060, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 17643, 3, 11152,
                                                                       11180, 15096, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 17688, 3, 11180,
                                                                       11208, 15132, ncols,
                                                                       gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 17733, 0, 3,
                                                                       17508, 14988, 17553,
                                                                       11264, 11348, 15168,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 17868, 0, 3,
                                                                       17553, 15024, 17598,
                                                                       11348, 11432, 15276,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 18003, 0, 3,
                                                                       17598, 15060, 17643,
                                                                       11432, 11516, 15384,
                                                                       ncols, gamma, p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 18138, 0, 3,
                                                                       17643, 15096, 17688,
                                                                       11516, 11600, 15492,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 18273, 0, 3,
                                                                       17733, 15168, 17868,
                                                                       11768, 11936, 15600,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 18543, 0, 3,
                                                                       17868, 15276, 18003,
                                                                       11936, 12104, 15816,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 18813, 0, 3,
                                                                       18003, 15384, 18138,
                                                                       12104, 12272, 16032,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 19083, 0, 3,
                                                                       18273, 15600, 18543,
                                                                       12608, 12888, 16248,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsl_three_center_electron_repulsion_0(buffer, 19533, 0, 3,
                                                                       18543, 15816, 18813,
                                                                       12888, 13168, 16608,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsl_three_center_electron_repulsion_0(buffer, 19983, 0, 3,
                                                                       19083, 16248, 19533,
                                                                       13728, 14148, 16968,
                                                                       ncols, gamma, p, q);

                    simdfunc::contract_primitives(buffer, 20658, 18273, 270, ncols);

                    simdfunc::contract_primitives(buffer, 21030, 19083, 450, ncols);

                    simdfunc::contract_primitives(buffer, 21650, 19983, 675, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 20928, 20658, 6, 1, nmax);

        simdtrf::transform_l_inner(buffer, 21480, 21030, 10, 1, nmax);

        simdtrf::transform_l_inner(buffer, 22325, 21650, 15, 1, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 22580, 20928, 21480, 17,
                                             nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 22886, 21480, 22325, 17,
                                             nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 23396, 22580, 22886, 17,
                                             nmax);

        simdtrf::transform_d_inner(buffer, 24008, 23396, 6, 17, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 24008, 85, nmax);
    }

    for (size_t m = 0; m < 425; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
