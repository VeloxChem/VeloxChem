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


#include "SimdThreeCenterElectronRepulsionRecGPH.hpp"

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
#include "SimdTransferGP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_gph_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_gph_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 11894, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 297 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 11894, 9752, 921, dimensions);

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 7, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10}, ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 18, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 21, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 24, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 27, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 8, 9,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 9, 10,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 10, 11,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 63, 0, 3, 11, 12,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 69, 0, 3, 12, 13,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 75, 0, 3, 13, 14,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 81, 0, 3, 14, 15,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 87, 0, 3, 15, 16,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 93, 0, 3, 18, 21,
                                                                       45, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 103, 0, 3, 21, 24,
                                                                       51, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 113, 0, 3, 24, 27,
                                                                       57, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 123, 0, 3, 27, 30,
                                                                       63, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 133, 0, 3, 30, 33,
                                                                       69, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 143, 0, 3, 33, 36,
                                                                       75, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 153, 0, 3, 36, 39,
                                                                       81, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 163, 0, 3, 45, 51,
                                                                       93, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 178, 0, 3, 51, 57,
                                                                       103, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 193, 0, 3, 57, 63,
                                                                       113, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 208, 0, 3, 63, 69,
                                                                       123, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 223, 0, 3, 69, 75,
                                                                       133, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 75, 81,
                                                                       143, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 253, 0, 3, 93,
                                                                       103, 163, 178, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 274, 0, 3, 103,
                                                                       113, 178, 193, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 295, 0, 3, 113,
                                                                       123, 193, 208, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 316, 0, 3, 123,
                                                                       133, 208, 223, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 337, 0, 3, 133,
                                                                       143, 223, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 358, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 361, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 364, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 367, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 370, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 373, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 376, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 379, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 382, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 385, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 388, 3, 8, 18,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 397, 3, 9, 21,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 406, 3, 10, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 415, 3, 11, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 424, 3, 12, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 433, 3, 13, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 442, 3, 14, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 451, 3, 15, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 460, 3, 16, 42,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 469, 3, 18, 45,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 487, 3, 21, 51,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 505, 3, 24, 57,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 523, 3, 27, 63,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 541, 3, 30, 69,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 559, 3, 33, 75,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 577, 3, 36, 81,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 595, 3, 39, 87,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 613, 3, 45, 93,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 643, 3, 51, 103,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 673, 3, 57, 113,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 703, 3, 63, 123,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 733, 3, 69, 133,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 763, 3, 75, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 793, 3, 81, 153,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 823, 3, 93, 163,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 868, 3, 103, 178,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 913, 3, 113, 193,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 958, 3, 123, 208,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1003, 3, 133, 223,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1048, 3, 143, 238,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1093, 3, 163, 253,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1156, 3, 178, 274,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1219, 3, 193, 295,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1282, 3, 208, 316,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1345, 3, 223, 337,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1408, 3, 8, 9,
                                                                       364, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1414, 3, 9, 10,
                                                                       367, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1420, 3, 10, 11,
                                                                       370, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1426, 3, 11, 12,
                                                                       373, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1432, 3, 12, 13,
                                                                       376, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1438, 3, 13, 14,
                                                                       379, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1444, 3, 14, 15,
                                                                       382, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1450, 3, 15, 16,
                                                                       385, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1456, 0, 3, 1408,
                                                                       364, 1414, 406, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1474, 0, 3, 1414,
                                                                       367, 1420, 415, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 1420,
                                                                       370, 1426, 424, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1510, 0, 3, 1426,
                                                                       373, 1432, 433, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1528, 0, 3, 1432,
                                                                       376, 1438, 442, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1546, 0, 3, 1438,
                                                                       379, 1444, 451, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1564, 0, 3, 1444,
                                                                       382, 1450, 460, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1582, 0, 3, 1456,
                                                                       406, 1474, 45, 51, 505,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1618, 0, 3, 1474,
                                                                       415, 1492, 51, 57, 523,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1654, 0, 3, 1492,
                                                                       424, 1510, 57, 63, 541,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1690, 0, 3, 1510,
                                                                       433, 1528, 63, 69, 559,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1726, 0, 3, 1528,
                                                                       442, 1546, 69, 75, 577,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1762, 0, 3, 1546,
                                                                       451, 1564, 75, 81, 595,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1798, 0, 3, 1582,
                                                                       505, 1618, 93, 103, 673,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1858, 0, 3, 1618,
                                                                       523, 1654, 103, 113, 703,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1918, 0, 3, 1654,
                                                                       541, 1690, 113, 123, 733,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1978, 0, 3, 1690,
                                                                       559, 1726, 123, 133, 763,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2038, 0, 3, 1726,
                                                                       577, 1762, 133, 143, 793,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2098, 0, 3, 1798,
                                                                       673, 1858, 163, 178, 913,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2188, 0, 3, 1858,
                                                                       703, 1918, 178, 193, 958,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2278, 0, 3, 1918,
                                                                       733, 1978, 193, 208, 1003,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2368, 0, 3, 1978,
                                                                       763, 2038, 208, 223, 1048,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 2458, 0, 3, 2098,
                                                                       913, 2188, 253, 274, 1219,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 2584, 0, 3, 2188,
                                                                       958, 2278, 274, 295, 1282,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 2710, 0, 3, 2278,
                                                                       1003, 2368, 295, 316,
                                                                       1345, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2836, 3, 358, 361,
                                                                       1408, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2846, 3, 361, 364,
                                                                       1414, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2856, 3, 364, 367,
                                                                       1420, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2866, 3, 367, 370,
                                                                       1426, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2876, 3, 370, 373,
                                                                       1432, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2886, 3, 373, 376,
                                                                       1438, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2896, 3, 376, 379,
                                                                       1444, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2906, 3, 379, 382,
                                                                       1450, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2916, 0, 3, 2836,
                                                                       1408, 2846, 388, 397,
                                                                       1456, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2946, 0, 3, 2846,
                                                                       1414, 2856, 397, 406,
                                                                       1474, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2976, 0, 3, 2856,
                                                                       1420, 2866, 406, 415,
                                                                       1492, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3006, 0, 3, 2866,
                                                                       1426, 2876, 415, 424,
                                                                       1510, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3036, 0, 3, 2876,
                                                                       1432, 2886, 424, 433,
                                                                       1528, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3066, 0, 3, 2886,
                                                                       1438, 2896, 433, 442,
                                                                       1546, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3096, 0, 3, 2896,
                                                                       1444, 2906, 442, 451,
                                                                       1564, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3126, 0, 3, 2916,
                                                                       1456, 2946, 469, 487,
                                                                       1582, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3186, 0, 3, 2946,
                                                                       1474, 2976, 487, 505,
                                                                       1618, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3246, 0, 3, 2976,
                                                                       1492, 3006, 505, 523,
                                                                       1654, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3306, 0, 3, 3006,
                                                                       1510, 3036, 523, 541,
                                                                       1690, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3366, 0, 3, 3036,
                                                                       1528, 3066, 541, 559,
                                                                       1726, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3426, 0, 3, 3066,
                                                                       1546, 3096, 559, 577,
                                                                       1762, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3486, 0, 3, 3126,
                                                                       1582, 3186, 613, 643,
                                                                       1798, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3586, 0, 3, 3186,
                                                                       1618, 3246, 643, 673,
                                                                       1858, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3686, 0, 3, 3246,
                                                                       1654, 3306, 673, 703,
                                                                       1918, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3786, 0, 3, 3306,
                                                                       1690, 3366, 703, 733,
                                                                       1978, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3886, 0, 3, 3366,
                                                                       1726, 3426, 733, 763,
                                                                       2038, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 3986, 0, 3, 3486,
                                                                       1798, 3586, 823, 868,
                                                                       2098, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4136, 0, 3, 3586,
                                                                       1858, 3686, 868, 913,
                                                                       2188, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4286, 0, 3, 3686,
                                                                       1918, 3786, 913, 958,
                                                                       2278, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4436, 0, 3, 3786,
                                                                       1978, 3886, 958, 1003,
                                                                       2368, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 4586, 0, 3, 3986,
                                                                       2098, 4136, 1093, 1156,
                                                                       2458, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 4796, 0, 3, 4136,
                                                                       2188, 4286, 1156, 1219,
                                                                       2584, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 5006, 0, 3, 4286,
                                                                       2278, 4436, 1219, 1282,
                                                                       2710, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5216, 3, 1408,
                                                                       1414, 2856, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5231, 3, 1414,
                                                                       1420, 2866, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5246, 3, 1420,
                                                                       1426, 2876, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5261, 3, 1426,
                                                                       1432, 2886, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5276, 3, 1432,
                                                                       1438, 2896, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 5291, 3, 1438,
                                                                       1444, 2906, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5306, 0, 3, 5216,
                                                                       2856, 5231, 1456, 1474,
                                                                       2976, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5351, 0, 3, 5231,
                                                                       2866, 5246, 1474, 1492,
                                                                       3006, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5396, 0, 3, 5246,
                                                                       2876, 5261, 1492, 1510,
                                                                       3036, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5441, 0, 3, 5261,
                                                                       2886, 5276, 1510, 1528,
                                                                       3066, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5486, 0, 3, 5276,
                                                                       2896, 5291, 1528, 1546,
                                                                       3096, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5531, 0, 3, 5306,
                                                                       2976, 5351, 1582, 1618,
                                                                       3246, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5621, 0, 3, 5351,
                                                                       3006, 5396, 1618, 1654,
                                                                       3306, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5711, 0, 3, 5396,
                                                                       3036, 5441, 1654, 1690,
                                                                       3366, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5801, 0, 3, 5441,
                                                                       3066, 5486, 1690, 1726,
                                                                       3426, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 5891, 0, 3, 5531,
                                                                       3246, 5621, 1798, 1858,
                                                                       3686, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6041, 0, 3, 5621,
                                                                       3306, 5711, 1858, 1918,
                                                                       3786, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 6191, 0, 3, 5711,
                                                                       3366, 5801, 1918, 1978,
                                                                       3886, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 6341, 0, 3, 5891,
                                                                       3686, 6041, 2098, 2188,
                                                                       4286, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 6566, 0, 3, 6041,
                                                                       3786, 6191, 2188, 2278,
                                                                       4436, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 6791, 0, 3, 6341,
                                                                       4286, 6566, 2458, 2584,
                                                                       5006, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7106, 3, 2836,
                                                                       2846, 5216, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7127, 3, 2846,
                                                                       2856, 5231, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7148, 3, 2856,
                                                                       2866, 5246, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7169, 3, 2866,
                                                                       2876, 5261, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7190, 3, 2876,
                                                                       2886, 5276, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 7211, 3, 2886,
                                                                       2896, 5291, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7232, 0, 3, 7106,
                                                                       5216, 7127, 2916, 2946,
                                                                       5306, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7295, 0, 3, 7127,
                                                                       5231, 7148, 2946, 2976,
                                                                       5351, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7358, 0, 3, 7148,
                                                                       5246, 7169, 2976, 3006,
                                                                       5396, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7421, 0, 3, 7169,
                                                                       5261, 7190, 3006, 3036,
                                                                       5441, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 7484, 0, 3, 7190,
                                                                       5276, 7211, 3036, 3066,
                                                                       5486, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 7547, 0, 3, 7232,
                                                                       5306, 7295, 3126, 3186,
                                                                       5531, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 7673, 0, 3, 7295,
                                                                       5351, 7358, 3186, 3246,
                                                                       5621, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 7799, 0, 3, 7358,
                                                                       5396, 7421, 3246, 3306,
                                                                       5711, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 7925, 0, 3, 7421,
                                                                       5441, 7484, 3306, 3366,
                                                                       5801, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 8051, 0, 3, 7547,
                                                                       5531, 7673, 3486, 3586,
                                                                       5891, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 8261, 0, 3, 7673,
                                                                       5621, 7799, 3586, 3686,
                                                                       6041, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 8471, 0, 3, 7799,
                                                                       5711, 7925, 3686, 3786,
                                                                       6191, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 8681, 0, 3, 8051,
                                                                       5891, 8261, 3986, 4136,
                                                                       6341, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 8996, 0, 3, 8261,
                                                                       6041, 8471, 4136, 4286,
                                                                       6566, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 9311, 0, 3, 8681,
                                                                       6341, 8996, 4586, 4796,
                                                                       6791, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 9752, 8681, 315, ncols);

                    simdfunc::contract_primitives(buffer, 10232, 9311, 441, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 10067, 9752, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 10673, 10232, 21, 1, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 10904, 10067, 10673, 11,
                                             nmax);

        simdtrf::transform_p_inner(buffer, 11399, 10904, 15, 11, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 11399, 33, nmax);
    }

    for (size_t m = 0; m < 297; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
