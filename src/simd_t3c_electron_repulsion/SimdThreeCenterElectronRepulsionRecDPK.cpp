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


#include "SimdThreeCenterElectronRepulsionRecDPK.hpp"

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
#include "SimdTransferDP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformK.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_dpk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_dpk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 8764, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 225 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 8764, 7408, 666, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 163, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 166, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 169, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 172, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 175, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 178, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 181, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 184, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 187, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 190, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 193, 3, 8, 18,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 202, 3, 9, 21,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 211, 3, 10, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 220, 3, 11, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 229, 3, 12, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 238, 3, 13, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 247, 3, 14, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 256, 3, 15, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 265, 3, 16, 42,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 274, 3, 18, 45,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 292, 3, 21, 51,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 310, 3, 24, 57,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 328, 3, 27, 63,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 346, 3, 30, 69,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 364, 3, 33, 75,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 382, 3, 36, 81,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 400, 3, 39, 87,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 418, 3, 45, 93,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 448, 3, 51, 103,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 478, 3, 57, 113,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 508, 3, 63, 123,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 538, 3, 69, 133,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 568, 3, 75, 143,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 598, 3, 81, 153,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 628, 3, 8, 9, 169,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 634, 3, 9, 10,
                                                                       172, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 640, 3, 10, 11,
                                                                       175, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 646, 3, 11, 12,
                                                                       178, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 652, 3, 12, 13,
                                                                       181, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 658, 3, 13, 14,
                                                                       184, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 664, 3, 14, 15,
                                                                       187, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 670, 3, 15, 16,
                                                                       190, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 676, 0, 3, 628,
                                                                       169, 634, 211, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 694, 0, 3, 634,
                                                                       172, 640, 220, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 712, 0, 3, 640,
                                                                       175, 646, 229, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 730, 0, 3, 646,
                                                                       178, 652, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 748, 0, 3, 652,
                                                                       181, 658, 247, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 766, 0, 3, 658,
                                                                       184, 664, 256, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 784, 0, 3, 664,
                                                                       187, 670, 265, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 802, 0, 3, 676,
                                                                       211, 694, 45, 51, 310,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 838, 0, 3, 694,
                                                                       220, 712, 51, 57, 328,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 874, 0, 3, 712,
                                                                       229, 730, 57, 63, 346,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 910, 0, 3, 730,
                                                                       238, 748, 63, 69, 364,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 946, 0, 3, 748,
                                                                       247, 766, 69, 75, 382,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 982, 0, 3, 766,
                                                                       256, 784, 75, 81, 400,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1018, 0, 3, 802,
                                                                       310, 838, 93, 103, 478,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1078, 0, 3, 838,
                                                                       328, 874, 103, 113, 508,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1138, 0, 3, 874,
                                                                       346, 910, 113, 123, 538,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1198, 0, 3, 910,
                                                                       364, 946, 123, 133, 568,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1258, 0, 3, 946,
                                                                       382, 982, 133, 143, 598,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1318, 3, 163, 166,
                                                                       628, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1328, 3, 166, 169,
                                                                       634, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1338, 3, 169, 172,
                                                                       640, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1348, 3, 172, 175,
                                                                       646, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1358, 3, 175, 178,
                                                                       652, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1368, 3, 178, 181,
                                                                       658, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1378, 3, 181, 184,
                                                                       664, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1388, 3, 184, 187,
                                                                       670, ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1398, 0, 3, 1318,
                                                                       628, 1328, 193, 202, 676,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1428, 0, 3, 1328,
                                                                       634, 1338, 202, 211, 694,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1458, 0, 3, 1338,
                                                                       640, 1348, 211, 220, 712,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1488, 0, 3, 1348,
                                                                       646, 1358, 220, 229, 730,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1518, 0, 3, 1358,
                                                                       652, 1368, 229, 238, 748,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 1368,
                                                                       658, 1378, 238, 247, 766,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1578, 0, 3, 1378,
                                                                       664, 1388, 247, 256, 784,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1608, 0, 3, 1398,
                                                                       676, 1428, 274, 292, 802,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1668, 0, 3, 1428,
                                                                       694, 1458, 292, 310, 838,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1728, 0, 3, 1458,
                                                                       712, 1488, 310, 328, 874,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1788, 0, 3, 1488,
                                                                       730, 1518, 328, 346, 910,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1848, 0, 3, 1518,
                                                                       748, 1548, 346, 364, 946,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1908, 0, 3, 1548,
                                                                       766, 1578, 364, 382, 982,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 1968, 0, 3, 1608,
                                                                       802, 1668, 418, 448, 1018,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 2068, 0, 3, 1668,
                                                                       838, 1728, 448, 478, 1078,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 2168, 0, 3, 1728,
                                                                       874, 1788, 478, 508, 1138,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 2268, 0, 3, 1788,
                                                                       910, 1848, 508, 538, 1198,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 2368, 0, 3, 1848,
                                                                       946, 1908, 538, 568, 1258,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2468, 3, 628, 634,
                                                                       1338, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2483, 3, 634, 640,
                                                                       1348, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2498, 3, 640, 646,
                                                                       1358, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2513, 3, 646, 652,
                                                                       1368, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2528, 3, 652, 658,
                                                                       1378, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2543, 3, 658, 664,
                                                                       1388, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 2558, 0, 3, 2468,
                                                                       1338, 2483, 676, 694,
                                                                       1458, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 2603, 0, 3, 2483,
                                                                       1348, 2498, 694, 712,
                                                                       1488, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 2498,
                                                                       1358, 2513, 712, 730,
                                                                       1518, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 2693, 0, 3, 2513,
                                                                       1368, 2528, 730, 748,
                                                                       1548, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 2738, 0, 3, 2528,
                                                                       1378, 2543, 748, 766,
                                                                       1578, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 2783, 0, 3, 2558,
                                                                       1458, 2603, 802, 838,
                                                                       1728, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 2873, 0, 3, 2603,
                                                                       1488, 2648, 838, 874,
                                                                       1788, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 2963, 0, 3, 2648,
                                                                       1518, 2693, 874, 910,
                                                                       1848, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 3053, 0, 3, 2693,
                                                                       1548, 2738, 910, 946,
                                                                       1908, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 3143, 0, 3, 2783,
                                                                       1728, 2873, 1018, 1078,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 3293, 0, 3, 2873,
                                                                       1788, 2963, 1078, 1138,
                                                                       2268, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 3443, 0, 3, 2963,
                                                                       1848, 3053, 1138, 1198,
                                                                       2368, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3593, 3, 1318,
                                                                       1328, 2468, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3614, 3, 1328,
                                                                       1338, 2483, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3635, 3, 1338,
                                                                       1348, 2498, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3656, 3, 1348,
                                                                       1358, 2513, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3677, 3, 1358,
                                                                       1368, 2528, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3698, 3, 1368,
                                                                       1378, 2543, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 3719, 0, 3, 3593,
                                                                       2468, 3614, 1398, 1428,
                                                                       2558, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 3782, 0, 3, 3614,
                                                                       2483, 3635, 1428, 1458,
                                                                       2603, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 3845, 0, 3, 3635,
                                                                       2498, 3656, 1458, 1488,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 3908, 0, 3, 3656,
                                                                       2513, 3677, 1488, 1518,
                                                                       2693, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 3971, 0, 3, 3677,
                                                                       2528, 3698, 1518, 1548,
                                                                       2738, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 4034, 0, 3, 3719,
                                                                       2558, 3782, 1608, 1668,
                                                                       2783, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 4160, 0, 3, 3782,
                                                                       2603, 3845, 1668, 1728,
                                                                       2873, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 4286, 0, 3, 3845,
                                                                       2648, 3908, 1728, 1788,
                                                                       2963, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 4412, 0, 3, 3908,
                                                                       2693, 3971, 1788, 1848,
                                                                       3053, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 4538, 0, 3, 4034,
                                                                       2783, 4160, 1968, 2068,
                                                                       3143, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 4748, 0, 3, 4160,
                                                                       2873, 4286, 2068, 2168,
                                                                       3293, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 4958, 0, 3, 4286,
                                                                       2963, 4412, 2168, 2268,
                                                                       3443, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5168, 3, 2468,
                                                                       2483, 3635, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5196, 3, 2483,
                                                                       2498, 3656, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5224, 3, 2498,
                                                                       2513, 3677, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5252, 3, 2513,
                                                                       2528, 3698, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 5280, 0, 3, 5168,
                                                                       3635, 5196, 2558, 2603,
                                                                       3845, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 5364, 0, 3, 5196,
                                                                       3656, 5224, 2603, 2648,
                                                                       3908, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 5448, 0, 3, 5224,
                                                                       3677, 5252, 2648, 2693,
                                                                       3971, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 5532, 0, 3, 5280,
                                                                       3845, 5364, 2783, 2873,
                                                                       4286, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 5700, 0, 3, 5364,
                                                                       3908, 5448, 2873, 2963,
                                                                       4412, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 5868, 0, 3, 5532,
                                                                       4286, 5700, 3143, 3293,
                                                                       4958, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 6148, 3, 3593,
                                                                       3614, 5168, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 6184, 3, 3614,
                                                                       3635, 5196, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 6220, 3, 3635,
                                                                       3656, 5224, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 6256, 3, 3656,
                                                                       3677, 5252, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 6292, 0, 3, 6148,
                                                                       5168, 6184, 3719, 3782,
                                                                       5280, ncols, gamma, p,
                                                                       q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 6400, 0, 3, 6184,
                                                                       5196, 6220, 3782, 3845,
                                                                       5364, ncols, gamma, p,
                                                                       q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 6508, 0, 3, 6220,
                                                                       5224, 6256, 3845, 3908,
                                                                       5448, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 6616, 0, 3, 6292,
                                                                       5280, 6400, 4034, 4160,
                                                                       5532, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 6832, 0, 3, 6400,
                                                                       5364, 6508, 4160, 4286,
                                                                       5700, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsk_three_center_electron_repulsion_0(buffer, 7048, 0, 3, 6616,
                                                                       5532, 6832, 4538, 4748,
                                                                       5868, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 7408, 6616, 216, ncols);

                    simdfunc::contract_primitives(buffer, 7714, 7048, 360, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 7624, 7408, 6, 1, nmax);

        simdtrf::transform_k_inner(buffer, 8074, 7714, 10, 1, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 8224, 7624, 8074, 15, nmax);

        simdtrf::transform_p_inner(buffer, 8494, 8224, 6, 15, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 8494, 45, nmax);
    }

    for (size_t m = 0; m < 225; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
