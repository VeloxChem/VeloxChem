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


#include "SimdThreeCenterElectronRepulsionRecDSL.hpp"

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
#include "SimdTransformD.hpp"
#include "SimdTransformL.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_dsl_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_dsl_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 5455, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 85 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 5455, 5083, 270, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 7, 3, 10,
                                                             ncols, fj, 6, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 19, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 22, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 8, 9,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 9, 10,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 61, 0, 3, 10, 11,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 67, 0, 3, 11, 12,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 73, 0, 3, 12, 13,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 79, 0, 3, 13, 14,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 85, 0, 3, 14, 15,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 91, 0, 3, 15, 16,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 97, 0, 3, 16, 17,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 103, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 106, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 109, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 112, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 115, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 118, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 121, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 124, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 127, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 130, 3, 10, 25,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 139, 3, 11, 28,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 148, 3, 12, 31,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 157, 3, 13, 34,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 166, 3, 14, 37,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 175, 3, 15, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 184, 3, 16, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 193, 3, 17, 46,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 202, 3, 25, 61,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 220, 3, 28, 67,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 238, 3, 31, 73,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 256, 3, 34, 79,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 274, 3, 37, 85,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 292, 3, 40, 91,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 310, 3, 43, 97,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 328, 3, 8, 9, 103,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 334, 3, 9, 10,
                                                                       106, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 340, 3, 10, 11,
                                                                       109, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 346, 3, 11, 12,
                                                                       112, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 352, 3, 12, 13,
                                                                       115, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 358, 3, 13, 14,
                                                                       118, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 364, 3, 14, 15,
                                                                       121, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 370, 3, 15, 16,
                                                                       124, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 376, 3, 16, 17,
                                                                       127, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 382, 0, 3, 328,
                                                                       103, 334, 130, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 400, 0, 3, 334,
                                                                       106, 340, 139, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 418, 0, 3, 340,
                                                                       109, 346, 148, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 436, 0, 3, 346,
                                                                       112, 352, 157, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 454, 0, 3, 352,
                                                                       115, 358, 166, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 472, 0, 3, 358,
                                                                       118, 364, 175, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 490, 0, 3, 364,
                                                                       121, 370, 184, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 508, 0, 3, 370,
                                                                       124, 376, 193, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 526, 0, 3, 382,
                                                                       130, 400, 49, 55, 202,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 562, 0, 3, 400,
                                                                       139, 418, 55, 61, 220,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 598, 0, 3, 418,
                                                                       148, 436, 61, 67, 238,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 634, 0, 3, 436,
                                                                       157, 454, 67, 73, 256,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 670, 0, 3, 454,
                                                                       166, 472, 73, 79, 274,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 706, 0, 3, 472,
                                                                       175, 490, 79, 85, 292,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 742, 0, 3, 490,
                                                                       184, 508, 85, 91, 310,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 778, 3, 103, 106,
                                                                       340, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 788, 3, 106, 109,
                                                                       346, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 798, 3, 109, 112,
                                                                       352, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 808, 3, 112, 115,
                                                                       358, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 818, 3, 115, 118,
                                                                       364, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 828, 3, 118, 121,
                                                                       370, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 838, 3, 121, 124,
                                                                       376, ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 848, 0, 3, 778,
                                                                       340, 788, 130, 139, 418,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 878, 0, 3, 788,
                                                                       346, 798, 139, 148, 436,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 908, 0, 3, 798,
                                                                       352, 808, 148, 157, 454,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 938, 0, 3, 808,
                                                                       358, 818, 157, 166, 472,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 968, 0, 3, 818,
                                                                       364, 828, 166, 175, 490,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 998, 0, 3, 828,
                                                                       370, 838, 175, 184, 508,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1028, 0, 3, 848,
                                                                       418, 878, 202, 220, 598,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1088, 0, 3, 878,
                                                                       436, 908, 220, 238, 634,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1148, 0, 3, 908,
                                                                       454, 938, 238, 256, 670,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1208, 0, 3, 938,
                                                                       472, 968, 256, 274, 706,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 968,
                                                                       490, 998, 274, 292, 742,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1328, 3, 328, 334,
                                                                       778, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1343, 3, 334, 340,
                                                                       788, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1358, 3, 340, 346,
                                                                       798, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1373, 3, 346, 352,
                                                                       808, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1388, 3, 352, 358,
                                                                       818, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1403, 3, 358, 364,
                                                                       828, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1418, 3, 364, 370,
                                                                       838, ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1433, 0, 3, 1328,
                                                                       778, 1343, 382, 400, 848,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1478, 0, 3, 1343,
                                                                       788, 1358, 400, 418, 878,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1523, 0, 3, 1358,
                                                                       798, 1373, 418, 436, 908,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1568, 0, 3, 1373,
                                                                       808, 1388, 436, 454, 938,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1613, 0, 3, 1388,
                                                                       818, 1403, 454, 472, 968,
                                                                       ncols, gamma, p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 1658, 0, 3, 1403,
                                                                       828, 1418, 472, 490, 998,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 1703, 0, 3, 1433,
                                                                       848, 1478, 526, 562, 1028,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 1793, 0, 3, 1478,
                                                                       878, 1523, 562, 598, 1088,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 1883, 0, 3, 1523,
                                                                       908, 1568, 598, 634, 1148,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 1973, 0, 3, 1568,
                                                                       938, 1613, 634, 670, 1208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 2063, 0, 3, 1613,
                                                                       968, 1658, 670, 706, 1268,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 2153, 3, 778, 788,
                                                                       1358, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 2174, 3, 788, 798,
                                                                       1373, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 2195, 3, 798, 808,
                                                                       1388, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 2216, 3, 808, 818,
                                                                       1403, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 2237, 3, 818, 828,
                                                                       1418, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 2258, 0, 3, 2153,
                                                                       1358, 2174, 848, 878,
                                                                       1523, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 2321, 0, 3, 2174,
                                                                       1373, 2195, 878, 908,
                                                                       1568, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 2384, 0, 3, 2195,
                                                                       1388, 2216, 908, 938,
                                                                       1613, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 2447, 0, 3, 2216,
                                                                       1403, 2237, 938, 968,
                                                                       1658, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 2510, 0, 3, 2258,
                                                                       1523, 2321, 1028, 1088,
                                                                       1883, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 2636, 0, 3, 2321,
                                                                       1568, 2384, 1088, 1148,
                                                                       1973, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 2762, 0, 3, 2384,
                                                                       1613, 2447, 1148, 1208,
                                                                       2063, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 2888, 3, 1328,
                                                                       1343, 2153, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 2916, 3, 1343,
                                                                       1358, 2174, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 2944, 3, 1358,
                                                                       1373, 2195, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 2972, 3, 1373,
                                                                       1388, 2216, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 3000, 3, 1388,
                                                                       1403, 2237, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 3028, 0, 3, 2888,
                                                                       2153, 2916, 1433, 1478,
                                                                       2258, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 3112, 0, 3, 2916,
                                                                       2174, 2944, 1478, 1523,
                                                                       2321, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 3196, 0, 3, 2944,
                                                                       2195, 2972, 1523, 1568,
                                                                       2384, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 3280, 0, 3, 2972,
                                                                       2216, 3000, 1568, 1613,
                                                                       2447, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 3364, 0, 3, 3028,
                                                                       2258, 3112, 1703, 1793,
                                                                       2510, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 3532, 0, 3, 3112,
                                                                       2321, 3196, 1793, 1883,
                                                                       2636, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 3700, 0, 3, 3196,
                                                                       2384, 3280, 1883, 1973,
                                                                       2762, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 3868, 3, 2153,
                                                                       2174, 2944, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 3904, 3, 2174,
                                                                       2195, 2972, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 3940, 3, 2195,
                                                                       2216, 3000, ncols, gamma,
                                                                       p, q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 3976, 0, 3, 3868,
                                                                       2944, 3904, 2258, 2321,
                                                                       3196, ncols, gamma, p,
                                                                       q);

                    compute_prim_psk_three_center_electron_repulsion_0(buffer, 4084, 0, 3, 3904,
                                                                       2972, 3940, 2321, 2384,
                                                                       3280, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsk_three_center_electron_repulsion_0(buffer, 4192, 0, 3, 3976,
                                                                       3196, 4084, 2510, 2636,
                                                                       3700, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 4408, 3, 2888,
                                                                       2916, 3868, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 4453, 3, 2916,
                                                                       2944, 3904, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssl_three_center_electron_repulsion_0(buffer, 4498, 3, 2944,
                                                                       2972, 3940, ncols, gamma,
                                                                       p, q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 4543, 0, 3, 4408,
                                                                       3868, 4453, 3028, 3112,
                                                                       3976, ncols, gamma, p,
                                                                       q);

                    compute_prim_psl_three_center_electron_repulsion_0(buffer, 4678, 0, 3, 4453,
                                                                       3904, 4498, 3112, 3196,
                                                                       4084, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsl_three_center_electron_repulsion_0(buffer, 4813, 0, 3, 4543,
                                                                       3976, 4678, 3364, 3532,
                                                                       4192, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 5083, 4813, 270, ncols);
                }
            }
        }

        simdtrf::transform_l_inner(buffer, 5353, 5083, 6, 1, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 5353, 17, nmax);
    }

    for (size_t m = 0; m < 85; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
