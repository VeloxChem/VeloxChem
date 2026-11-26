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


#include "SimdThreeCenterElectronRepulsionGeom010RecDSI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_dsi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_dsi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 10656, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 195 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 10656, 10074, 504, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 9, 6, 9, ncols,
                                                             fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 20, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 56, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 59, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 62, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 65, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 68, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 71, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 74, 0, 6, 10, 11,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 83, 0, 6, 11, 12,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 92, 0, 6, 12, 13,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 101, 0, 6, 13, 14,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 110, 0, 6, 14, 15,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 119, 0, 6, 15, 16,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 128, 0, 6, 16, 17,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 137, 0, 6, 17, 18,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 146, 0, 6, 10, 11,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 152, 0, 6, 11, 12,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 158, 0, 6, 12, 13,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 164, 0, 6, 13, 14,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 170, 0, 6, 14, 15,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 176, 0, 6, 15, 16,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 182, 0, 6, 16, 17,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 188, 0, 6, 17, 18,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 194, 0, 3, 6, 47,
                                                                       50, 74, 83, 146, 152,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 212, 0, 3, 6, 50,
                                                                       53, 83, 92, 152, 158,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 230, 0, 3, 6, 53,
                                                                       56, 92, 101, 158, 164,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 248, 0, 3, 6, 56,
                                                                       59, 101, 110, 164, 170,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 266, 0, 3, 6, 59,
                                                                       62, 110, 119, 170, 176,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 284, 0, 3, 6, 62,
                                                                       65, 119, 128, 176, 182,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 302, 0, 3, 6, 65,
                                                                       68, 128, 137, 182, 188,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 320, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 323, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 326, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 329, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 332, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 335, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 338, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 341, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 344, 6, 12, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 353, 6, 13, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 362, 6, 14, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 371, 6, 15, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 380, 6, 16, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 389, 6, 17, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 398, 6, 18, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 407, 6, 12, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 416, 6, 13, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 425, 6, 14, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 434, 6, 15, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 443, 6, 16, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 452, 6, 17, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 461, 6, 18, 71,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 470, 6, 26, 53,
                                                                       92, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 497, 6, 29, 56,
                                                                       101, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 524, 6, 32, 59,
                                                                       110, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 551, 6, 35, 62,
                                                                       119, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 578, 6, 38, 65,
                                                                       128, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 605, 6, 41, 68,
                                                                       137, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 632, 6, 53, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 650, 6, 56, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 668, 6, 59, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 686, 6, 62, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 704, 6, 65, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 722, 6, 68, 188,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 740, 0, 6, 470,
                                                                       92, 497, 158, 230, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 794, 0, 6, 497,
                                                                       101, 524, 164, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 848, 0, 6, 524,
                                                                       110, 551, 170, 266, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 902, 0, 6, 551,
                                                                       119, 578, 176, 284, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 956, 0, 6, 578,
                                                                       128, 605, 182, 302, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1010, 6, 10, 11,
                                                                       320, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1016, 6, 11, 12,
                                                                       323, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1022, 6, 12, 13,
                                                                       326, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1028, 6, 13, 14,
                                                                       329, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1034, 6, 14, 15,
                                                                       332, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1040, 6, 15, 16,
                                                                       335, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1046, 6, 16, 17,
                                                                       338, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1052, 6, 17, 18,
                                                                       341, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1058, 3, 6, 1010,
                                                                       320, 1016, 344, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1076, 3, 6, 1016,
                                                                       323, 1022, 353, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1094, 3, 6, 1022,
                                                                       326, 1028, 362, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1112, 3, 6, 1028,
                                                                       329, 1034, 371, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1130, 3, 6, 1034,
                                                                       332, 1040, 380, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1148, 3, 6, 1040,
                                                                       335, 1046, 389, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1166, 3, 6, 1046,
                                                                       338, 1052, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1184, 0, 6, 1010,
                                                                       320, 1016, 407, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1202, 0, 6, 1016,
                                                                       323, 1022, 416, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1220, 0, 6, 1022,
                                                                       326, 1028, 425, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1238, 0, 6, 1028,
                                                                       329, 1034, 434, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1256, 0, 6, 1034,
                                                                       332, 1040, 443, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1274, 0, 6, 1040,
                                                                       335, 1046, 452, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1292, 0, 6, 1046,
                                                                       338, 1052, 461, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1310, 0, 3, 6,
                                                                       1058, 344, 1076, 1184,
                                                                       407, 1202, 74, 83, 470,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1364, 0, 3, 6,
                                                                       1076, 353, 1094, 1202,
                                                                       416, 1220, 83, 92, 497,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1418, 0, 3, 6,
                                                                       1094, 362, 1112, 1220,
                                                                       425, 1238, 92, 101, 524,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1472, 0, 3, 6,
                                                                       1112, 371, 1130, 1238,
                                                                       434, 1256, 101, 110, 551,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1526, 0, 3, 6,
                                                                       1130, 380, 1148, 1256,
                                                                       443, 1274, 110, 119, 578,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1580, 0, 3, 6,
                                                                       1148, 389, 1166, 1274,
                                                                       452, 1292, 119, 128, 605,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1634, 0, 6, 1184,
                                                                       407, 1202, 146, 152, 632,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1670, 0, 6, 1202,
                                                                       416, 1220, 152, 158, 650,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1706, 0, 6, 1220,
                                                                       425, 1238, 158, 164, 668,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1742, 0, 6, 1238,
                                                                       434, 1256, 164, 170, 686,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1778, 0, 6, 1256,
                                                                       443, 1274, 170, 176, 704,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1814, 0, 6, 1274,
                                                                       452, 1292, 176, 182, 722,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 1850, 0, 3, 6,
                                                                       1310, 470, 1364, 1634,
                                                                       632, 1670, 194, 212, 740,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 1958, 0, 3, 6,
                                                                       1364, 497, 1418, 1670,
                                                                       650, 1706, 212, 230, 794,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 2066, 0, 3, 6,
                                                                       1418, 524, 1472, 1706,
                                                                       668, 1742, 230, 248, 848,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 2174, 0, 3, 6,
                                                                       1472, 551, 1526, 1742,
                                                                       686, 1778, 248, 266, 902,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 2282, 0, 3, 6,
                                                                       1526, 578, 1580, 1778,
                                                                       704, 1814, 266, 284, 956,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2390, 6, 320, 323,
                                                                       1022, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2400, 6, 323, 326,
                                                                       1028, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2410, 6, 326, 329,
                                                                       1034, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2420, 6, 329, 332,
                                                                       1040, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2430, 6, 332, 335,
                                                                       1046, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 2440, 6, 335, 338,
                                                                       1052, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2450, 3, 6, 2390,
                                                                       1022, 2400, 1094, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2480, 3, 6, 2400,
                                                                       1028, 2410, 1112, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2510, 3, 6, 2410,
                                                                       1034, 2420, 1130, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2540, 3, 6, 2420,
                                                                       1040, 2430, 1148, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2570, 3, 6, 2430,
                                                                       1046, 2440, 1166, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2600, 0, 6, 2390,
                                                                       1022, 2400, 407, 416,
                                                                       1220, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2630, 0, 6, 2400,
                                                                       1028, 2410, 416, 425,
                                                                       1238, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2660, 0, 6, 2410,
                                                                       1034, 2420, 425, 434,
                                                                       1256, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2690, 0, 6, 2420,
                                                                       1040, 2430, 434, 443,
                                                                       1274, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2720, 0, 6, 2430,
                                                                       1046, 2440, 443, 452,
                                                                       1292, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2750, 0, 3, 6,
                                                                       2450, 1094, 2480, 2600,
                                                                       1220, 2630, 470, 497,
                                                                       1418, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2840, 0, 3, 6,
                                                                       2480, 1112, 2510, 2630,
                                                                       1238, 2660, 497, 524,
                                                                       1472, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2930, 0, 3, 6,
                                                                       2510, 1130, 2540, 2660,
                                                                       1256, 2690, 524, 551,
                                                                       1526, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 3020, 0, 3, 6,
                                                                       2540, 1148, 2570, 2690,
                                                                       1274, 2720, 551, 578,
                                                                       1580, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3110, 0, 6, 2600,
                                                                       1220, 2630, 632, 650,
                                                                       1706, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3170, 0, 6, 2630,
                                                                       1238, 2660, 650, 668,
                                                                       1742, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3230, 0, 6, 2660,
                                                                       1256, 2690, 668, 686,
                                                                       1778, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3290, 0, 6, 2690,
                                                                       1274, 2720, 686, 704,
                                                                       1814, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 3350, 0, 3, 6,
                                                                       2750, 1418, 2840, 3110,
                                                                       1706, 3170, 740, 794,
                                                                       2066, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 3530, 0, 3, 6,
                                                                       2840, 1472, 2930, 3170,
                                                                       1742, 3230, 794, 848,
                                                                       2174, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 3710, 0, 3, 6,
                                                                       2930, 1526, 3020, 3230,
                                                                       1778, 3290, 848, 902,
                                                                       2282, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3890, 6, 1010,
                                                                       1016, 2390, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3905, 6, 1016,
                                                                       1022, 2400, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3920, 6, 1022,
                                                                       1028, 2410, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3935, 6, 1028,
                                                                       1034, 2420, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3950, 6, 1034,
                                                                       1040, 2430, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3965, 6, 1040,
                                                                       1046, 2440, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3980, 3, 6, 3890,
                                                                       2390, 3905, 1058, 1076,
                                                                       2450, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4025, 3, 6, 3905,
                                                                       2400, 3920, 1076, 1094,
                                                                       2480, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4070, 3, 6, 3920,
                                                                       2410, 3935, 1094, 1112,
                                                                       2510, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4115, 3, 6, 3935,
                                                                       2420, 3950, 1112, 1130,
                                                                       2540, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4160, 3, 6, 3950,
                                                                       2430, 3965, 1130, 1148,
                                                                       2570, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4205, 0, 6, 3890,
                                                                       2390, 3905, 1184, 1202,
                                                                       2600, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4250, 0, 6, 3905,
                                                                       2400, 3920, 1202, 1220,
                                                                       2630, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4295, 0, 6, 3920,
                                                                       2410, 3935, 1220, 1238,
                                                                       2660, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4340, 0, 6, 3935,
                                                                       2420, 3950, 1238, 1256,
                                                                       2690, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 4385, 0, 6, 3950,
                                                                       2430, 3965, 1256, 1274,
                                                                       2720, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 4430, 0, 3, 6,
                                                                       3980, 2450, 4025, 4205,
                                                                       2600, 4250, 1310, 1364,
                                                                       2750, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 4565, 0, 3, 6,
                                                                       4025, 2480, 4070, 4250,
                                                                       2630, 4295, 1364, 1418,
                                                                       2840, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 4700, 0, 3, 6,
                                                                       4070, 2510, 4115, 4295,
                                                                       2660, 4340, 1418, 1472,
                                                                       2930, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 4835, 0, 3, 6,
                                                                       4115, 2540, 4160, 4340,
                                                                       2690, 4385, 1472, 1526,
                                                                       3020, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 4970, 0, 6, 4205,
                                                                       2600, 4250, 1634, 1670,
                                                                       3110, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5060, 0, 6, 4250,
                                                                       2630, 4295, 1670, 1706,
                                                                       3170, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5150, 0, 6, 4295,
                                                                       2660, 4340, 1706, 1742,
                                                                       3230, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 5240, 0, 6, 4340,
                                                                       2690, 4385, 1742, 1778,
                                                                       3290, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 5330, 0, 3, 6,
                                                                       4430, 2750, 4565, 4970,
                                                                       3110, 5060, 1850, 1958,
                                                                       3350, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 5600, 0, 3, 6,
                                                                       4565, 2840, 4700, 5060,
                                                                       3170, 5150, 1958, 2066,
                                                                       3530, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 5870, 0, 3, 6,
                                                                       4700, 2930, 4835, 5150,
                                                                       3230, 5240, 2066, 2174,
                                                                       3710, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6140, 6, 2390,
                                                                       2400, 3920, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6161, 6, 2400,
                                                                       2410, 3935, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6182, 6, 2410,
                                                                       2420, 3950, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 6203, 6, 2420,
                                                                       2430, 3965, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 6224, 3, 6, 6140,
                                                                       3920, 6161, 2450, 2480,
                                                                       4070, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 6287, 3, 6, 6161,
                                                                       3935, 6182, 2480, 2510,
                                                                       4115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 6350, 3, 6, 6182,
                                                                       3950, 6203, 2510, 2540,
                                                                       4160, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 6413, 0, 6, 6140,
                                                                       3920, 6161, 2600, 2630,
                                                                       4295, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 6476, 0, 6, 6161,
                                                                       3935, 6182, 2630, 2660,
                                                                       4340, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 6539, 0, 6, 6182,
                                                                       3950, 6203, 2660, 2690,
                                                                       4385, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 6602, 0, 3, 6,
                                                                       6224, 4070, 6287, 6413,
                                                                       4295, 6476, 2750, 2840,
                                                                       4700, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 6791, 0, 3, 6,
                                                                       6287, 4115, 6350, 6476,
                                                                       4340, 6539, 2840, 2930,
                                                                       4835, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 6980, 0, 6, 6413,
                                                                       4295, 6476, 3110, 3170,
                                                                       5150, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 7106, 0, 6, 6476,
                                                                       4340, 6539, 3170, 3230,
                                                                       5240, ncols, gamma, p,
                                                                       q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 7232, 0, 3, 6,
                                                                       6602, 4700, 6791, 6980,
                                                                       5150, 7106, 3350, 3530,
                                                                       5870, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 7610, 6, 3890,
                                                                       3905, 6140, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 7638, 6, 3905,
                                                                       3920, 6161, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 7666, 6, 3920,
                                                                       3935, 6182, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 7694, 6, 3935,
                                                                       3950, 6203, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 7722, 3, 6, 7610,
                                                                       6140, 7638, 3980, 4025,
                                                                       6224, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 7806, 3, 6, 7638,
                                                                       6161, 7666, 4025, 4070,
                                                                       6287, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 7890, 3, 6, 7666,
                                                                       6182, 7694, 4070, 4115,
                                                                       6350, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 7974, 0, 6, 7610,
                                                                       6140, 7638, 4205, 4250,
                                                                       6413, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 8058, 0, 6, 7638,
                                                                       6161, 7666, 4250, 4295,
                                                                       6476, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 8142, 0, 6, 7666,
                                                                       6182, 7694, 4295, 4340,
                                                                       6539, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 8226, 0, 3, 6,
                                                                       7722, 6224, 7806, 7974,
                                                                       6413, 8058, 4430, 4565,
                                                                       6602, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 8478, 0, 3, 6,
                                                                       7806, 6287, 7890, 8058,
                                                                       6476, 8142, 4565, 4700,
                                                                       6791, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 8730, 0, 6, 7974,
                                                                       6413, 8058, 4970, 5060,
                                                                       6980, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 8898, 0, 6, 8058,
                                                                       6476, 8142, 5060, 5150,
                                                                       7106, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 9066, 0, 3, 6,
                                                                       8226, 6602, 8478, 8730,
                                                                       6980, 8898, 5330, 5600,
                                                                       7232, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_s_x(buffer, 9570, 9066, 6, 28, ncols, beta);

                    simdgeo::geom_s_y(buffer, 9738, 9066, 6, 28, ncols, beta);

                    simdgeo::geom_s_z(buffer, 9906, 9066, 6, 28, ncols, beta);

                    simdfunc::contract_primitives(buffer, 10074, 9570, 504, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 10578, 10074, 6, 1, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 10578, 13, nmax);

        simdtrf::transform_i_inner(buffer, 10578, 10242, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 65 * nvalues + n * npairs, nvalues, buffer, 10578,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 10578, 10410, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 130 * nvalues + n * npairs, nvalues, buffer, 10578,
                                   13, nmax);
    }

    for (size_t m = 0; m < 195; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
