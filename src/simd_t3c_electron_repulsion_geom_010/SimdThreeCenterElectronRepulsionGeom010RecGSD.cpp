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


#include "SimdThreeCenterElectronRepulsionGeom010RecGSD.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdGeometryS1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_gsd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_gsd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    // NOTE: a derivative screens with the integral's own bound, on purpose. It
    // is not a bound on the derivative -- that is larger by roughly 2 alpha R,
    // the relation reaching one shell higher -- and tightening it here would be
    // the wrong repair. A screened Fock build defines an energy in which the
    // dropped pairs contribute exactly zero, and the derivative of that energy
    // is the derivative screened the same way; a tighter bound would add forces
    // from pairs the energy never counted. One threshold controls both errors,
    // so tightening it in the Fock build tightens the gradient with it.

    const auto dimensions = simdfunc::make_column_dimensions(
        a_function, b_function, c_function, npairs, coordinates,
        screenfunc::three_center_electron_repulsion_primitive_bound,
        threshold / static_cast<double>(nprims));

    const auto nmax = simdfunc::prepare_buffer(buffer, 4055, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 135 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 4055, 3710, 270, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto beta = b_exps[j];

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fa = -b_exps[j] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pa(buffer, coordinates, 0, nmax, fa);

                simdfunc::compute_pb(buffer, coordinates, 3, nmax, fb);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 6, n, nmax, fc);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 9, 6, 7, ncols,
                                                             fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 18, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 60, 0, 6, 10, 11,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 69, 0, 6, 11, 12,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 78, 0, 6, 12, 13,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 87, 0, 6, 13, 14,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 96, 0, 6, 14, 15,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 105, 0, 6, 15, 16,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 114, 0, 6, 10, 11,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 120, 0, 6, 11, 12,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 126, 0, 6, 12, 13,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 132, 0, 6, 13, 14,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 138, 0, 6, 14, 15,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 144, 0, 6, 15, 16,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 150, 0, 3, 6, 39,
                                                                       42, 60, 69, 114, 120,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 168, 0, 3, 6, 42,
                                                                       45, 69, 78, 120, 126,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 186, 0, 3, 6, 45,
                                                                       48, 78, 87, 126, 132,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 204, 0, 3, 6, 48,
                                                                       51, 87, 96, 132, 138,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 222, 0, 3, 6, 51,
                                                                       54, 96, 105, 138, 144,
                                                                       ncols, gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 240, 0, 6, 39, 42,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 250, 0, 6, 42, 45,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 260, 0, 6, 45, 48,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 270, 0, 6, 48, 51,
                                                                       132, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 280, 0, 6, 51, 54,
                                                                       138, 144, ncols, gamma, p,
                                                                       q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 290, 0, 3, 6, 114,
                                                                       120, 150, 168, 240, 250,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 320, 0, 3, 6, 120,
                                                                       126, 168, 186, 250, 260,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 350, 0, 3, 6, 126,
                                                                       132, 186, 204, 260, 270,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 380, 0, 3, 6, 132,
                                                                       138, 204, 222, 270, 280,
                                                                       ncols, gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 410, 0, 6, 114,
                                                                       120, 240, 250, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 425, 0, 6, 120,
                                                                       126, 250, 260, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 440, 0, 6, 126,
                                                                       132, 260, 270, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 455, 0, 6, 132,
                                                                       138, 270, 280, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 470, 0, 3, 6, 150,
                                                                       168, 240, 250, 290, 320,
                                                                       410, 425, ncols, gamma, p,
                                                                       q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 515, 0, 3, 6, 168,
                                                                       186, 250, 260, 320, 350,
                                                                       425, 440, ncols, gamma, p,
                                                                       q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 560, 0, 3, 6, 186,
                                                                       204, 260, 270, 350, 380,
                                                                       440, 455, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 605, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 608, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 611, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 614, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 617, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 620, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 623, 6, 12, 24,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 632, 6, 13, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 641, 6, 14, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 650, 6, 15, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 659, 6, 16, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 668, 6, 12, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 677, 6, 13, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 686, 6, 14, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 695, 6, 15, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 704, 6, 16, 57,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 713, 6, 24, 45,
                                                                       78, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 740, 6, 27, 48,
                                                                       87, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 767, 6, 30, 51,
                                                                       96, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 794, 6, 33, 54,
                                                                       105, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 821, 6, 45, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 839, 6, 48, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 857, 6, 51, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 875, 6, 54, 144,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 893, 0, 6, 713,
                                                                       78, 740, 126, 186, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 947, 0, 6, 740,
                                                                       87, 767, 132, 204, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1001, 0, 6, 767,
                                                                       96, 794, 138, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1055, 6, 126, 260,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1085, 6, 132, 270,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1115, 6, 138, 280,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1145, 0, 6, 893,
                                                                       186, 947, 260, 350, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1235, 0, 6, 947,
                                                                       204, 1001, 270, 380,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1325, 6, 260, 440,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1370, 6, 270, 455,
                                                                       ncols, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 1415, 0, 6, 1145,
                                                                       350, 1235, 440, 560,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1550, 6, 10, 11,
                                                                       605, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1556, 6, 11, 12,
                                                                       608, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1562, 6, 12, 13,
                                                                       611, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1568, 6, 13, 14,
                                                                       614, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1574, 6, 14, 15,
                                                                       617, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1580, 6, 15, 16,
                                                                       620, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1586, 3, 6, 1550,
                                                                       605, 1556, 623, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1604, 3, 6, 1556,
                                                                       608, 1562, 632, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1622, 3, 6, 1562,
                                                                       611, 1568, 641, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1640, 3, 6, 1568,
                                                                       614, 1574, 650, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1658, 3, 6, 1574,
                                                                       617, 1580, 659, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1676, 0, 6, 1550,
                                                                       605, 1556, 668, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1694, 0, 6, 1556,
                                                                       608, 1562, 677, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1712, 0, 6, 1562,
                                                                       611, 1568, 686, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1730, 0, 6, 1568,
                                                                       614, 1574, 695, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1748, 0, 6, 1574,
                                                                       617, 1580, 704, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1766, 0, 3, 6,
                                                                       1586, 623, 1604, 1676,
                                                                       668, 1694, 60, 69, 713,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1820, 0, 3, 6,
                                                                       1604, 632, 1622, 1694,
                                                                       677, 1712, 69, 78, 740,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1874, 0, 3, 6,
                                                                       1622, 641, 1640, 1712,
                                                                       686, 1730, 78, 87, 767,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1928, 0, 3, 6,
                                                                       1640, 650, 1658, 1730,
                                                                       695, 1748, 87, 96, 794,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1982, 0, 6, 1676,
                                                                       668, 1694, 114, 120, 821,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2018, 0, 6, 1694,
                                                                       677, 1712, 120, 126, 839,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2054, 0, 6, 1712,
                                                                       686, 1730, 126, 132, 857,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2090, 0, 6, 1730,
                                                                       695, 1748, 132, 138, 875,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 2126, 0, 3, 6,
                                                                       1766, 713, 1820, 1982,
                                                                       821, 2018, 150, 168, 893,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 2234, 0, 3, 6,
                                                                       1820, 740, 1874, 2018,
                                                                       839, 2054, 168, 186, 947,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 2342, 0, 3, 6,
                                                                       1874, 767, 1928, 2054,
                                                                       857, 2090, 186, 204, 1001,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2450, 0, 6, 1982,
                                                                       821, 2018, 240, 250, 1055,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2510, 0, 6, 2018,
                                                                       839, 2054, 250, 260, 1085,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2570, 0, 6, 2054,
                                                                       857, 2090, 260, 270, 1115,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 2630, 0, 3, 6,
                                                                       1766, 1820, 2126, 893,
                                                                       2234, 2450, 1055, 2510,
                                                                       290, 320, 1145, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 2810, 0, 3, 6,
                                                                       1820, 1874, 2234, 947,
                                                                       2342, 2510, 1085, 2570,
                                                                       320, 350, 1235, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2990, 0, 6, 2450,
                                                                       1055, 2510, 410, 425,
                                                                       1325, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3080, 0, 6, 2510,
                                                                       1085, 2570, 425, 440,
                                                                       1370, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 3170, 0, 3, 6,
                                                                       2126, 2234, 2630, 1145,
                                                                       2810, 2990, 1325, 3080,
                                                                       470, 515, 1415, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 3440, 3170, 15, 6, ncols, beta);

                    simdgeo::geom_s_y(buffer, 3530, 3170, 15, 6, ncols, beta);

                    simdgeo::geom_s_z(buffer, 3620, 3170, 15, 6, ncols, beta);

                    simdfunc::contract_primitives(buffer, 3710, 3440, 270, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 3980, 3710, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 3980, 5, nmax);

        simdtrf::transform_d_inner(buffer, 3980, 3800, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 45 * nvalues + n * npairs, nvalues, buffer, 3980, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 3980, 3890, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 90 * nvalues + n * npairs, nvalues, buffer, 3980, 5,
                                   nmax);
    }

    for (size_t m = 0; m < 135; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
