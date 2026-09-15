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


#include "SimdThreeCenterElectronRepulsionRecSDK.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSK.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformK.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_sdk_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_sdk_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 3696, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 75 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 3696, 3390, 216, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto mu = a_exps[i] * b_exps[j] / p;

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pb(buffer, coordinates, 0, nmax, fb);

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
                                                        5, 6, 7, 8, 9}, ncols, fj, 6, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 17, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 20, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 41, 0, 3, 8, 9,
                                                                       17, 20, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 47, 0, 3, 9, 10,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 53, 0, 3, 10, 11,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 59, 0, 3, 11, 12,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 65, 0, 3, 12, 13,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 71, 0, 3, 13, 14,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 77, 0, 3, 14, 15,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 83, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 86, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 89, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 92, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 95, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 98, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 101, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 104, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 107, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 110, 3, 10, 23,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 119, 3, 11, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 128, 3, 12, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 137, 3, 13, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 146, 3, 14, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 155, 3, 15, 38,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 164, 3, 17, 41,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 182, 3, 20, 47,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 200, 3, 23, 53,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 218, 3, 26, 59,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 236, 3, 29, 65,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 254, 3, 32, 71,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 272, 3, 35, 77,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 290, 3, 8, 9, 89,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 296, 3, 9, 10, 92,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 302, 3, 10, 11,
                                                                       95, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 308, 3, 11, 12,
                                                                       98, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 314, 3, 12, 13,
                                                                       101, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 320, 3, 13, 14,
                                                                       104, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 326, 3, 14, 15,
                                                                       107, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 332, 0, 3, 290,
                                                                       89, 296, 110, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 350, 0, 3, 296,
                                                                       92, 302, 119, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 368, 0, 3, 302,
                                                                       95, 308, 128, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 386, 0, 3, 308,
                                                                       98, 314, 137, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 404, 0, 3, 314,
                                                                       101, 320, 146, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 422, 0, 3, 320,
                                                                       104, 326, 155, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 440, 0, 3, 332,
                                                                       110, 350, 41, 47, 200,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 476, 0, 3, 350,
                                                                       119, 368, 47, 53, 218,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 512, 0, 3, 368,
                                                                       128, 386, 53, 59, 236,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 548, 0, 3, 386,
                                                                       137, 404, 59, 65, 254,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 584, 0, 3, 404,
                                                                       146, 422, 65, 71, 272,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 620, 3, 83, 86,
                                                                       290, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 630, 3, 86, 89,
                                                                       296, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 640, 3, 89, 92,
                                                                       302, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 650, 3, 92, 95,
                                                                       308, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 660, 3, 95, 98,
                                                                       314, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 670, 3, 98, 101,
                                                                       320, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 680, 3, 101, 104,
                                                                       326, ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 690, 0, 3, 620,
                                                                       290, 630, 332, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 720, 0, 3, 630,
                                                                       296, 640, 350, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 750, 0, 3, 640,
                                                                       302, 650, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 780, 0, 3, 650,
                                                                       308, 660, 386, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 810, 0, 3, 660,
                                                                       314, 670, 404, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 840, 0, 3, 670,
                                                                       320, 680, 422, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 870, 0, 3, 690,
                                                                       332, 720, 164, 182, 440,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 930, 0, 3, 720,
                                                                       350, 750, 182, 200, 476,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 990, 0, 3, 750,
                                                                       368, 780, 200, 218, 512,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1050, 0, 3, 780,
                                                                       386, 810, 218, 236, 548,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 1110, 0, 3, 810,
                                                                       404, 840, 236, 254, 584,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1170, 3, 290, 296,
                                                                       640, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1185, 3, 296, 302,
                                                                       650, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1200, 3, 302, 308,
                                                                       660, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1215, 3, 308, 314,
                                                                       670, ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 1230, 3, 314, 320,
                                                                       680, ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 1245, 0, 3, 1170,
                                                                       640, 1185, 332, 350, 750,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 1290, 0, 3, 1185,
                                                                       650, 1200, 350, 368, 780,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 1335, 0, 3, 1200,
                                                                       660, 1215, 368, 386, 810,
                                                                       ncols, gamma, p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 1215,
                                                                       670, 1230, 386, 404, 840,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 1425, 0, 3, 1245,
                                                                       750, 1290, 440, 476, 990,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 1515, 0, 3, 1290,
                                                                       780, 1335, 476, 512, 1050,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 1605, 0, 3, 1335,
                                                                       810, 1380, 512, 548, 1110,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 1695, 3, 620, 630,
                                                                       1170, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 1716, 3, 630, 640,
                                                                       1185, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 1737, 3, 640, 650,
                                                                       1200, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 1758, 3, 650, 660,
                                                                       1215, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 1779, 3, 660, 670,
                                                                       1230, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 1800, 0, 3, 1695,
                                                                       1170, 1716, 690, 720,
                                                                       1245, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 1863, 0, 3, 1716,
                                                                       1185, 1737, 720, 750,
                                                                       1290, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 1926, 0, 3, 1737,
                                                                       1200, 1758, 750, 780,
                                                                       1335, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 1989, 0, 3, 1758,
                                                                       1215, 1779, 780, 810,
                                                                       1380, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 2052, 0, 3, 1800,
                                                                       1245, 1863, 870, 930,
                                                                       1425, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 2178, 0, 3, 1863,
                                                                       1290, 1926, 930, 990,
                                                                       1515, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 2304, 0, 3, 1926,
                                                                       1335, 1989, 990, 1050,
                                                                       1605, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 2430, 3, 1170,
                                                                       1185, 1737, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 2458, 3, 1185,
                                                                       1200, 1758, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 2486, 3, 1200,
                                                                       1215, 1779, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 2514, 0, 3, 2430,
                                                                       1737, 2458, 1245, 1290,
                                                                       1926, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 2598, 0, 3, 2458,
                                                                       1758, 2486, 1290, 1335,
                                                                       1989, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 2682, 0, 3, 2514,
                                                                       1926, 2598, 1425, 1515,
                                                                       2304, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 2850, 3, 1695,
                                                                       1716, 2430, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 2886, 3, 1716,
                                                                       1737, 2458, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssk_three_center_electron_repulsion_0(buffer, 2922, 3, 1737,
                                                                       1758, 2486, ncols, gamma,
                                                                       p, q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 2958, 0, 3, 2850,
                                                                       2430, 2886, 1800, 1863,
                                                                       2514, ncols, gamma, p,
                                                                       q);

                    compute_prim_spk_three_center_electron_repulsion_0(buffer, 3066, 0, 3, 2886,
                                                                       2458, 2922, 1863, 1926,
                                                                       2598, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdk_three_center_electron_repulsion_0(buffer, 3174, 0, 3, 2958,
                                                                       2514, 3066, 2052, 2178,
                                                                       2682, ncols, gamma, p,
                                                                       q);

                    simdfunc::contract_primitives(buffer, 3390, 3174, 216, ncols);
                }
            }
        }

        simdtrf::transform_k_inner(buffer, 3606, 3390, 6, 1, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 3606, 15, nmax);
    }

    for (size_t m = 0; m < 75; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
