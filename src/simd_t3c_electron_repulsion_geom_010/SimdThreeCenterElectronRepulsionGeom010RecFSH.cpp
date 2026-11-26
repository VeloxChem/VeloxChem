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


#include "SimdThreeCenterElectronRepulsionGeom010RecFSH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_fsh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_fsh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 14162, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 231 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 14162, 13422, 630, dimensions);

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 9, 6, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 19, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 67, 0, 6, 10, 11,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 76, 0, 6, 11, 12,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 85, 0, 6, 12, 13,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 94, 0, 6, 13, 14,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 103, 0, 6, 14, 15,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 112, 0, 6, 15, 16,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 121, 0, 6, 16, 17,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 130, 0, 6, 10, 11,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 136, 0, 6, 11, 12,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 142, 0, 6, 12, 13,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 148, 0, 6, 13, 14,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 154, 0, 6, 14, 15,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 160, 0, 6, 15, 16,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 166, 0, 6, 16, 17,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 172, 0, 3, 6, 43,
                                                                       46, 67, 76, 130, 136,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 190, 0, 3, 6, 46,
                                                                       49, 76, 85, 136, 142,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 208, 0, 3, 6, 49,
                                                                       52, 85, 94, 142, 148,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 226, 0, 3, 6, 52,
                                                                       55, 94, 103, 148, 154,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 244, 0, 3, 6, 55,
                                                                       58, 103, 112, 154, 160,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 262, 0, 3, 6, 58,
                                                                       61, 112, 121, 160, 166,
                                                                       ncols, gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 280, 0, 6, 43, 46,
                                                                       130, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 290, 0, 6, 46, 49,
                                                                       136, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 300, 0, 6, 49, 52,
                                                                       142, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 310, 0, 6, 52, 55,
                                                                       148, 154, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 320, 0, 6, 55, 58,
                                                                       154, 160, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 330, 0, 6, 58, 61,
                                                                       160, 166, ncols, gamma, p,
                                                                       q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 340, 0, 3, 6, 130,
                                                                       136, 172, 190, 280, 290,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 370, 0, 3, 6, 136,
                                                                       142, 190, 208, 290, 300,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 400, 0, 3, 6, 142,
                                                                       148, 208, 226, 300, 310,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 430, 0, 3, 6, 148,
                                                                       154, 226, 244, 310, 320,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 460, 0, 3, 6, 154,
                                                                       160, 244, 262, 320, 330,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 490, 6, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 493, 6, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 496, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 499, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 502, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 505, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 508, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 511, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 514, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 517, 6, 12, 25,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 526, 6, 13, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 535, 6, 14, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 544, 6, 15, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 553, 6, 16, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 562, 6, 17, 40,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 571, 6, 10, 43,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 580, 6, 11, 46,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 589, 6, 12, 49,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 598, 6, 13, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 607, 6, 14, 55,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 616, 6, 15, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 625, 6, 16, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 634, 6, 17, 64,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 643, 6, 19, 43,
                                                                       67, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 670, 6, 22, 46,
                                                                       76, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 697, 6, 25, 49,
                                                                       85, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 724, 6, 28, 52,
                                                                       94, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 751, 6, 31, 55,
                                                                       103, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 778, 6, 34, 58,
                                                                       112, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 805, 6, 37, 61,
                                                                       121, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 832, 6, 43, 130,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 850, 6, 46, 136,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 868, 6, 49, 142,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 886, 6, 52, 148,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 904, 6, 55, 154,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 922, 6, 58, 160,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 940, 6, 61, 166,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 958, 0, 6, 643,
                                                                       67, 670, 130, 172, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1012, 0, 6, 670,
                                                                       76, 697, 136, 190, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1066, 0, 6, 697,
                                                                       85, 724, 142, 208, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1120, 0, 6, 724,
                                                                       94, 751, 148, 226, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1174, 0, 6, 751,
                                                                       103, 778, 154, 244, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1228, 0, 6, 778,
                                                                       112, 805, 160, 262, ncols,
                                                                       gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1282, 6, 130, 280,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1312, 6, 136, 290,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1342, 6, 142, 300,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1372, 6, 148, 310,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1402, 6, 154, 320,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1432, 6, 160, 330,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1462, 0, 6, 958,
                                                                       172, 1012, 280, 340,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1552, 0, 6, 1012,
                                                                       190, 1066, 290, 370,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1642, 0, 6, 1066,
                                                                       208, 1120, 300, 400,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1732, 0, 6, 1120,
                                                                       226, 1174, 310, 430,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1822, 0, 6, 1174,
                                                                       244, 1228, 320, 460,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1912, 6, 10, 11,
                                                                       496, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1918, 6, 11, 12,
                                                                       499, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1924, 6, 12, 13,
                                                                       502, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1930, 6, 13, 14,
                                                                       505, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1936, 6, 14, 15,
                                                                       508, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1942, 6, 15, 16,
                                                                       511, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1948, 6, 16, 17,
                                                                       514, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1954, 3, 6, 1912,
                                                                       496, 1918, 517, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1972, 3, 6, 1918,
                                                                       499, 1924, 526, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1990, 3, 6, 1924,
                                                                       502, 1930, 535, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2008, 3, 6, 1930,
                                                                       505, 1936, 544, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2026, 3, 6, 1936,
                                                                       508, 1942, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2044, 3, 6, 1942,
                                                                       511, 1948, 562, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2062, 0, 6, 1912,
                                                                       496, 1918, 589, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2080, 0, 6, 1918,
                                                                       499, 1924, 598, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2098, 0, 6, 1924,
                                                                       502, 1930, 607, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2116, 0, 6, 1930,
                                                                       505, 1936, 616, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2134, 0, 6, 1936,
                                                                       508, 1942, 625, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2152, 0, 6, 1942,
                                                                       511, 1948, 634, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2170, 0, 3, 6,
                                                                       1954, 517, 1972, 2062,
                                                                       589, 2080, 67, 76, 697,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2224, 0, 3, 6,
                                                                       1972, 526, 1990, 2080,
                                                                       598, 2098, 76, 85, 724,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2278, 0, 3, 6,
                                                                       1990, 535, 2008, 2098,
                                                                       607, 2116, 85, 94, 751,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2332, 0, 3, 6,
                                                                       2008, 544, 2026, 2116,
                                                                       616, 2134, 94, 103, 778,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2386, 0, 3, 6,
                                                                       2026, 553, 2044, 2134,
                                                                       625, 2152, 103, 112, 805,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2440, 0, 6, 2062,
                                                                       589, 2080, 130, 136, 868,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2476, 0, 6, 2080,
                                                                       598, 2098, 136, 142, 886,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2512, 0, 6, 2098,
                                                                       607, 2116, 142, 148, 904,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2548, 0, 6, 2116,
                                                                       616, 2134, 148, 154, 922,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2584, 0, 6, 2134,
                                                                       625, 2152, 154, 160, 940,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 2620, 0, 3, 6,
                                                                       2170, 697, 2224, 2440,
                                                                       868, 2476, 172, 190, 1066,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 2728, 0, 3, 6,
                                                                       2224, 724, 2278, 2476,
                                                                       886, 2512, 190, 208, 1120,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 2836, 0, 3, 6,
                                                                       2278, 751, 2332, 2512,
                                                                       904, 2548, 208, 226, 1174,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 2944, 0, 3, 6,
                                                                       2332, 778, 2386, 2548,
                                                                       922, 2584, 226, 244, 1228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3052, 0, 6, 2440,
                                                                       868, 2476, 280, 290, 1342,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3112, 0, 6, 2476,
                                                                       886, 2512, 290, 300, 1372,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3172, 0, 6, 2512,
                                                                       904, 2548, 300, 310, 1402,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3232, 0, 6, 2548,
                                                                       922, 2584, 310, 320, 1432,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 3292, 0, 3, 6,
                                                                       2170, 2224, 2620, 1066,
                                                                       2728, 3052, 1342, 3112,
                                                                       340, 370, 1642, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 3472, 0, 3, 6,
                                                                       2224, 2278, 2728, 1120,
                                                                       2836, 3112, 1372, 3172,
                                                                       370, 400, 1732, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 3652, 0, 3, 6,
                                                                       2278, 2332, 2836, 1174,
                                                                       2944, 3172, 1402, 3232,
                                                                       400, 430, 1822, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3832, 6, 490, 493,
                                                                       1912, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3842, 6, 493, 496,
                                                                       1918, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3852, 6, 496, 499,
                                                                       1924, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3862, 6, 499, 502,
                                                                       1930, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3872, 6, 502, 505,
                                                                       1936, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3882, 6, 505, 508,
                                                                       1942, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3892, 6, 508, 511,
                                                                       1948, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3902, 3, 6, 3832,
                                                                       1912, 3842, 1954, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3932, 3, 6, 3842,
                                                                       1918, 3852, 1972, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3962, 3, 6, 3852,
                                                                       1924, 3862, 1990, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3992, 3, 6, 3862,
                                                                       1930, 3872, 2008, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4022, 3, 6, 3872,
                                                                       1936, 3882, 2026, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4052, 3, 6, 3882,
                                                                       1942, 3892, 2044, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4082, 0, 6, 3832,
                                                                       1912, 3842, 571, 580,
                                                                       2062, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4112, 0, 6, 3842,
                                                                       1918, 3852, 580, 589,
                                                                       2080, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4142, 0, 6, 3852,
                                                                       1924, 3862, 589, 598,
                                                                       2098, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4172, 0, 6, 3862,
                                                                       1930, 3872, 598, 607,
                                                                       2116, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4202, 0, 6, 3872,
                                                                       1936, 3882, 607, 616,
                                                                       2134, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4232, 0, 6, 3882,
                                                                       1942, 3892, 616, 625,
                                                                       2152, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4262, 0, 3, 6,
                                                                       3902, 1954, 3932, 4082,
                                                                       2062, 4112, 643, 670,
                                                                       2170, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4352, 0, 3, 6,
                                                                       3932, 1972, 3962, 4112,
                                                                       2080, 4142, 670, 697,
                                                                       2224, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4442, 0, 3, 6,
                                                                       3962, 1990, 3992, 4142,
                                                                       2098, 4172, 697, 724,
                                                                       2278, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4532, 0, 3, 6,
                                                                       3992, 2008, 4022, 4172,
                                                                       2116, 4202, 724, 751,
                                                                       2332, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4622, 0, 3, 6,
                                                                       4022, 2026, 4052, 4202,
                                                                       2134, 4232, 751, 778,
                                                                       2386, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4712, 0, 6, 4082,
                                                                       2062, 4112, 832, 850,
                                                                       2440, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4772, 0, 6, 4112,
                                                                       2080, 4142, 850, 868,
                                                                       2476, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4832, 0, 6, 4142,
                                                                       2098, 4172, 868, 886,
                                                                       2512, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4892, 0, 6, 4172,
                                                                       2116, 4202, 886, 904,
                                                                       2548, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4952, 0, 6, 4202,
                                                                       2134, 4232, 904, 922,
                                                                       2584, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 5012, 0, 3, 6,
                                                                       4262, 2170, 4352, 4712,
                                                                       2440, 4772, 958, 1012,
                                                                       2620, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 5192, 0, 3, 6,
                                                                       4352, 2224, 4442, 4772,
                                                                       2476, 4832, 1012, 1066,
                                                                       2728, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 5372, 0, 3, 6,
                                                                       4442, 2278, 4532, 4832,
                                                                       2512, 4892, 1066, 1120,
                                                                       2836, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 5552, 0, 3, 6,
                                                                       4532, 2332, 4622, 4892,
                                                                       2548, 4952, 1120, 1174,
                                                                       2944, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5732, 0, 6, 4712,
                                                                       2440, 4772, 1282, 1312,
                                                                       3052, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5832, 0, 6, 4772,
                                                                       2476, 4832, 1312, 1342,
                                                                       3112, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5932, 0, 6, 4832,
                                                                       2512, 4892, 1342, 1372,
                                                                       3172, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6032, 0, 6, 4892,
                                                                       2548, 4952, 1372, 1402,
                                                                       3232, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 6132, 0, 3, 6,
                                                                       4262, 4352, 5012, 2620,
                                                                       5192, 5732, 3052, 5832,
                                                                       1462, 1552, 3292, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 6432, 0, 3, 6,
                                                                       4352, 4442, 5192, 2728,
                                                                       5372, 5832, 3112, 5932,
                                                                       1552, 1642, 3472, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 6732, 0, 3, 6,
                                                                       4442, 4532, 5372, 2836,
                                                                       5552, 5932, 3172, 6032,
                                                                       1642, 1732, 3652, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7032, 6, 1912,
                                                                       1918, 3852, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7047, 6, 1918,
                                                                       1924, 3862, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7062, 6, 1924,
                                                                       1930, 3872, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7077, 6, 1930,
                                                                       1936, 3882, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7092, 6, 1936,
                                                                       1942, 3892, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7107, 3, 6, 7032,
                                                                       3852, 7047, 1954, 1972,
                                                                       3962, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7152, 3, 6, 7047,
                                                                       3862, 7062, 1972, 1990,
                                                                       3992, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7197, 3, 6, 7062,
                                                                       3872, 7077, 1990, 2008,
                                                                       4022, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7242, 3, 6, 7077,
                                                                       3882, 7092, 2008, 2026,
                                                                       4052, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7287, 0, 6, 7032,
                                                                       3852, 7047, 2062, 2080,
                                                                       4142, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7332, 0, 6, 7047,
                                                                       3862, 7062, 2080, 2098,
                                                                       4172, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7377, 0, 6, 7062,
                                                                       3872, 7077, 2098, 2116,
                                                                       4202, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7422, 0, 6, 7077,
                                                                       3882, 7092, 2116, 2134,
                                                                       4232, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 7467, 0, 3, 6,
                                                                       7107, 3962, 7152, 7287,
                                                                       4142, 7332, 2170, 2224,
                                                                       4442, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 7602, 0, 3, 6,
                                                                       7152, 3992, 7197, 7332,
                                                                       4172, 7377, 2224, 2278,
                                                                       4532, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 7737, 0, 3, 6,
                                                                       7197, 4022, 7242, 7377,
                                                                       4202, 7422, 2278, 2332,
                                                                       4622, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 7872, 0, 6, 7287,
                                                                       4142, 7332, 2440, 2476,
                                                                       4832, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 7962, 0, 6, 7332,
                                                                       4172, 7377, 2476, 2512,
                                                                       4892, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8052, 0, 6, 7377,
                                                                       4202, 7422, 2512, 2548,
                                                                       4952, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 8142, 0, 3, 6,
                                                                       7467, 4442, 7602, 7872,
                                                                       4832, 7962, 2620, 2728,
                                                                       5372, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 8412, 0, 3, 6,
                                                                       7602, 4532, 7737, 7962,
                                                                       4892, 8052, 2728, 2836,
                                                                       5552, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 8682, 0, 6, 7872,
                                                                       4832, 7962, 3052, 3112,
                                                                       5932, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 8832, 0, 6, 7962,
                                                                       4892, 8052, 3112, 3172,
                                                                       6032, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 8982, 0, 3, 6,
                                                                       7467, 7602, 8142, 5372,
                                                                       8412, 8682, 5932, 8832,
                                                                       3292, 3472, 6732, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9432, 6, 3832,
                                                                       3842, 7032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9453, 6, 3842,
                                                                       3852, 7047, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9474, 6, 3852,
                                                                       3862, 7062, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9495, 6, 3862,
                                                                       3872, 7077, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9516, 6, 3872,
                                                                       3882, 7092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9537, 3, 6, 9432,
                                                                       7032, 9453, 3902, 3932,
                                                                       7107, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9600, 3, 6, 9453,
                                                                       7047, 9474, 3932, 3962,
                                                                       7152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9663, 3, 6, 9474,
                                                                       7062, 9495, 3962, 3992,
                                                                       7197, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9726, 3, 6, 9495,
                                                                       7077, 9516, 3992, 4022,
                                                                       7242, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9789, 0, 6, 9432,
                                                                       7032, 9453, 4082, 4112,
                                                                       7287, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9852, 0, 6, 9453,
                                                                       7047, 9474, 4112, 4142,
                                                                       7332, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9915, 0, 6, 9474,
                                                                       7062, 9495, 4142, 4172,
                                                                       7377, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9978, 0, 6, 9495,
                                                                       7077, 9516, 4172, 4202,
                                                                       7422, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 10041, 0, 3, 6,
                                                                       9537, 7107, 9600, 9789,
                                                                       7287, 9852, 4262, 4352,
                                                                       7467, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 10230, 0, 3, 6,
                                                                       9600, 7152, 9663, 9852,
                                                                       7332, 9915, 4352, 4442,
                                                                       7602, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 10419, 0, 3, 6,
                                                                       9663, 7197, 9726, 9915,
                                                                       7377, 9978, 4442, 4532,
                                                                       7737, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 10608, 0, 6, 9789,
                                                                       7287, 9852, 4712, 4772,
                                                                       7872, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 10734, 0, 6, 9852,
                                                                       7332, 9915, 4772, 4832,
                                                                       7962, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 10860, 0, 6, 9915,
                                                                       7377, 9978, 4832, 4892,
                                                                       8052, ncols, gamma, p,
                                                                       q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 10986, 0, 3, 6,
                                                                       10041, 7467, 10230, 10608,
                                                                       7872, 10734, 5012, 5192,
                                                                       8142, ncols, gamma, p,
                                                                       q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 11364, 0, 3, 6,
                                                                       10230, 7602, 10419, 10734,
                                                                       7962, 10860, 5192, 5372,
                                                                       8412, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 11742, 0, 6,
                                                                       10608, 7872, 10734, 5732,
                                                                       5832, 8682, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 11952, 0, 6,
                                                                       10734, 7962, 10860, 5832,
                                                                       5932, 8832, ncols, gamma,
                                                                       p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 12162, 0, 3, 6,
                                                                       10041, 10230, 10986, 8142,
                                                                       11364, 11742, 8682, 11952,
                                                                       6132, 6432, 8982, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 12792, 12162, 10, 21, ncols, beta);

                    simdgeo::geom_s_y(buffer, 13002, 12162, 10, 21, ncols, beta);

                    simdgeo::geom_s_z(buffer, 13212, 12162, 10, 21, ncols, beta);

                    simdfunc::contract_primitives(buffer, 13422, 12792, 630, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 14052, 13422, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 14052, 11, nmax);

        simdtrf::transform_h_inner(buffer, 14052, 13632, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 77 * nvalues + n * npairs, nvalues, buffer, 14052,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 14052, 13842, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 154 * nvalues + n * npairs, nvalues, buffer, 14052,
                                   11, nmax);
    }

    for (size_t m = 0; m < 231; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
