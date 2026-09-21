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


#include "SimdThreeCenterElectronRepulsionGeom010RecFSI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
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
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_fsi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_fsi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 22620, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 273 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 22620, 21650, 840, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 9, 6, 10,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 3, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 60, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 63, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 66, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 69, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 72, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 75, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 78, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 81, 0, 6, 10, 11,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 90, 0, 6, 11, 12,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 99, 0, 6, 12, 13,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 108, 0, 6, 13, 14,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 117, 0, 6, 14, 15,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 126, 0, 6, 15, 16,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 135, 0, 6, 16, 17,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 144, 0, 6, 17, 18,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 153, 0, 6, 18, 19,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 162, 0, 6, 10, 11,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 168, 0, 6, 11, 12,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 174, 0, 6, 12, 13,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 180, 0, 6, 13, 14,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 186, 0, 6, 14, 15,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 192, 0, 6, 15, 16,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 198, 0, 6, 16, 17,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 204, 0, 6, 17, 18,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 210, 0, 6, 18, 19,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 216, 0, 3, 6, 51,
                                                                       54, 81, 90, 162, 168,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 234, 0, 3, 6, 54,
                                                                       57, 90, 99, 168, 174,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 252, 0, 3, 6, 57,
                                                                       60, 99, 108, 174, 180,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 270, 0, 3, 6, 60,
                                                                       63, 108, 117, 180, 186,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 288, 0, 3, 6, 63,
                                                                       66, 117, 126, 186, 192,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 306, 0, 3, 6, 66,
                                                                       69, 126, 135, 192, 198,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 324, 0, 3, 6, 69,
                                                                       72, 135, 144, 198, 204,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 342, 0, 3, 6, 72,
                                                                       75, 144, 153, 204, 210,
                                                                       ncols, gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 360, 0, 6, 51, 54,
                                                                       162, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 370, 0, 6, 54, 57,
                                                                       168, 174, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 380, 0, 6, 57, 60,
                                                                       174, 180, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 390, 0, 6, 60, 63,
                                                                       180, 186, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 400, 0, 6, 63, 66,
                                                                       186, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 410, 0, 6, 66, 69,
                                                                       192, 198, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 420, 0, 6, 69, 72,
                                                                       198, 204, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 430, 0, 6, 72, 75,
                                                                       204, 210, ncols, gamma, p,
                                                                       q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 440, 0, 3, 6, 162,
                                                                       168, 216, 234, 360, 370,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 470, 0, 3, 6, 168,
                                                                       174, 234, 252, 370, 380,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 500, 0, 3, 6, 174,
                                                                       180, 252, 270, 380, 390,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 530, 0, 3, 6, 180,
                                                                       186, 270, 288, 390, 400,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 560, 0, 3, 6, 186,
                                                                       192, 288, 306, 400, 410,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 590, 0, 3, 6, 192,
                                                                       198, 306, 324, 410, 420,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 620, 0, 3, 6, 198,
                                                                       204, 324, 342, 420, 430,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 650, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 653, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 656, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 659, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 662, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 665, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 668, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 671, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 674, 6, 20, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 677, 6, 12, 27,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 686, 6, 13, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 695, 6, 14, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 704, 6, 15, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 713, 6, 16, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 722, 6, 17, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 731, 6, 18, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 740, 6, 19, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 749, 6, 12, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 758, 6, 13, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 767, 6, 14, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 776, 6, 15, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 785, 6, 16, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 794, 6, 17, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 803, 6, 18, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 812, 6, 19, 78,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 821, 6, 27, 57,
                                                                       99, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 848, 6, 30, 60,
                                                                       108, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 875, 6, 33, 63,
                                                                       117, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 902, 6, 36, 66,
                                                                       126, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 929, 6, 39, 69,
                                                                       135, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 956, 6, 42, 72,
                                                                       144, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 983, 6, 45, 75,
                                                                       153, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1010, 6, 57, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1028, 6, 60, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1046, 6, 63, 186,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1064, 6, 66, 192,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1082, 6, 69, 198,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1100, 6, 72, 204,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1118, 6, 75, 210,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1136, 0, 6, 821,
                                                                       99, 848, 174, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1190, 0, 6, 848,
                                                                       108, 875, 180, 270, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1244, 0, 6, 875,
                                                                       117, 902, 186, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1298, 0, 6, 902,
                                                                       126, 929, 192, 306, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1352, 0, 6, 929,
                                                                       135, 956, 198, 324, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1406, 0, 6, 956,
                                                                       144, 983, 204, 342, ncols,
                                                                       gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1460, 6, 174, 380,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1490, 6, 180, 390,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1520, 6, 186, 400,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1550, 6, 192, 410,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1580, 6, 198, 420,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1610, 6, 204, 430,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1640, 0, 6, 1136,
                                                                       252, 1190, 380, 500,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1730, 0, 6, 1190,
                                                                       270, 1244, 390, 530,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1820, 0, 6, 1244,
                                                                       288, 1298, 400, 560,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1910, 0, 6, 1298,
                                                                       306, 1352, 410, 590,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2000, 0, 6, 1352,
                                                                       324, 1406, 420, 620,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2090, 6, 10, 11,
                                                                       650, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2096, 6, 11, 12,
                                                                       653, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2102, 6, 12, 13,
                                                                       656, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2108, 6, 13, 14,
                                                                       659, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2114, 6, 14, 15,
                                                                       662, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2120, 6, 15, 16,
                                                                       665, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2126, 6, 16, 17,
                                                                       668, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2132, 6, 17, 18,
                                                                       671, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2138, 6, 18, 19,
                                                                       674, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2144, 3, 6, 2090,
                                                                       650, 2096, 677, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2162, 3, 6, 2096,
                                                                       653, 2102, 686, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2180, 3, 6, 2102,
                                                                       656, 2108, 695, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2198, 3, 6, 2108,
                                                                       659, 2114, 704, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2216, 3, 6, 2114,
                                                                       662, 2120, 713, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2234, 3, 6, 2120,
                                                                       665, 2126, 722, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2252, 3, 6, 2126,
                                                                       668, 2132, 731, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2270, 3, 6, 2132,
                                                                       671, 2138, 740, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2288, 0, 6, 2090,
                                                                       650, 2096, 749, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2306, 0, 6, 2096,
                                                                       653, 2102, 758, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2324, 0, 6, 2102,
                                                                       656, 2108, 767, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2342, 0, 6, 2108,
                                                                       659, 2114, 776, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2360, 0, 6, 2114,
                                                                       662, 2120, 785, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2378, 0, 6, 2120,
                                                                       665, 2126, 794, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2396, 0, 6, 2126,
                                                                       668, 2132, 803, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2414, 0, 6, 2132,
                                                                       671, 2138, 812, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2432, 0, 3, 6,
                                                                       2144, 677, 2162, 2288,
                                                                       749, 2306, 81, 90, 821,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2486, 0, 3, 6,
                                                                       2162, 686, 2180, 2306,
                                                                       758, 2324, 90, 99, 848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2540, 0, 3, 6,
                                                                       2180, 695, 2198, 2324,
                                                                       767, 2342, 99, 108, 875,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2594, 0, 3, 6,
                                                                       2198, 704, 2216, 2342,
                                                                       776, 2360, 108, 117, 902,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2648, 0, 3, 6,
                                                                       2216, 713, 2234, 2360,
                                                                       785, 2378, 117, 126, 929,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2702, 0, 3, 6,
                                                                       2234, 722, 2252, 2378,
                                                                       794, 2396, 126, 135, 956,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2756, 0, 3, 6,
                                                                       2252, 731, 2270, 2396,
                                                                       803, 2414, 135, 144, 983,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2810, 0, 6, 2288,
                                                                       749, 2306, 162, 168, 1010,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2846, 0, 6, 2306,
                                                                       758, 2324, 168, 174, 1028,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2882, 0, 6, 2324,
                                                                       767, 2342, 174, 180, 1046,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2918, 0, 6, 2342,
                                                                       776, 2360, 180, 186, 1064,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2954, 0, 6, 2360,
                                                                       785, 2378, 186, 192, 1082,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2990, 0, 6, 2378,
                                                                       794, 2396, 192, 198, 1100,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3026, 0, 6, 2396,
                                                                       803, 2414, 198, 204, 1118,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3062, 0, 3, 6,
                                                                       2432, 821, 2486, 2810,
                                                                       1010, 2846, 216, 234,
                                                                       1136, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3170, 0, 3, 6,
                                                                       2486, 848, 2540, 2846,
                                                                       1028, 2882, 234, 252,
                                                                       1190, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3278, 0, 3, 6,
                                                                       2540, 875, 2594, 2882,
                                                                       1046, 2918, 252, 270,
                                                                       1244, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3386, 0, 3, 6,
                                                                       2594, 902, 2648, 2918,
                                                                       1064, 2954, 270, 288,
                                                                       1298, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3494, 0, 3, 6,
                                                                       2648, 929, 2702, 2954,
                                                                       1082, 2990, 288, 306,
                                                                       1352, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3602, 0, 3, 6,
                                                                       2702, 956, 2756, 2990,
                                                                       1100, 3026, 306, 324,
                                                                       1406, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3710, 0, 6, 2810,
                                                                       1010, 2846, 360, 370,
                                                                       1460, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3770, 0, 6, 2846,
                                                                       1028, 2882, 370, 380,
                                                                       1490, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3830, 0, 6, 2882,
                                                                       1046, 2918, 380, 390,
                                                                       1520, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3890, 0, 6, 2918,
                                                                       1064, 2954, 390, 400,
                                                                       1550, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3950, 0, 6, 2954,
                                                                       1082, 2990, 400, 410,
                                                                       1580, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4010, 0, 6, 2990,
                                                                       1100, 3026, 410, 420,
                                                                       1610, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 4070, 0, 3, 6,
                                                                       2432, 2486, 3062, 1136,
                                                                       3170, 3710, 1460, 3770,
                                                                       440, 470, 1640, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 4250, 0, 3, 6,
                                                                       2486, 2540, 3170, 1190,
                                                                       3278, 3770, 1490, 3830,
                                                                       470, 500, 1730, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 4430, 0, 3, 6,
                                                                       2540, 2594, 3278, 1244,
                                                                       3386, 3830, 1520, 3890,
                                                                       500, 530, 1820, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 4610, 0, 3, 6,
                                                                       2594, 2648, 3386, 1298,
                                                                       3494, 3890, 1550, 3950,
                                                                       530, 560, 1910, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 4790, 0, 3, 6,
                                                                       2648, 2702, 3494, 1352,
                                                                       3602, 3950, 1580, 4010,
                                                                       560, 590, 2000, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4970, 6, 650, 653,
                                                                       2102, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4980, 6, 653, 656,
                                                                       2108, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4990, 6, 656, 659,
                                                                       2114, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5000, 6, 659, 662,
                                                                       2120, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5010, 6, 662, 665,
                                                                       2126, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5020, 6, 665, 668,
                                                                       2132, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 5030, 6, 668, 671,
                                                                       2138, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5040, 3, 6, 4970,
                                                                       2102, 4980, 2180, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5070, 3, 6, 4980,
                                                                       2108, 4990, 2198, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5100, 3, 6, 4990,
                                                                       2114, 5000, 2216, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5130, 3, 6, 5000,
                                                                       2120, 5010, 2234, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5160, 3, 6, 5010,
                                                                       2126, 5020, 2252, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5190, 3, 6, 5020,
                                                                       2132, 5030, 2270, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5220, 0, 6, 4970,
                                                                       2102, 4980, 749, 758,
                                                                       2324, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5250, 0, 6, 4980,
                                                                       2108, 4990, 758, 767,
                                                                       2342, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5280, 0, 6, 4990,
                                                                       2114, 5000, 767, 776,
                                                                       2360, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5310, 0, 6, 5000,
                                                                       2120, 5010, 776, 785,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5340, 0, 6, 5010,
                                                                       2126, 5020, 785, 794,
                                                                       2396, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5370, 0, 6, 5020,
                                                                       2132, 5030, 794, 803,
                                                                       2414, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5400, 0, 3, 6,
                                                                       5040, 2180, 5070, 5220,
                                                                       2324, 5250, 821, 848,
                                                                       2540, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5490, 0, 3, 6,
                                                                       5070, 2198, 5100, 5250,
                                                                       2342, 5280, 848, 875,
                                                                       2594, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5580, 0, 3, 6,
                                                                       5100, 2216, 5130, 5280,
                                                                       2360, 5310, 875, 902,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5670, 0, 3, 6,
                                                                       5130, 2234, 5160, 5310,
                                                                       2378, 5340, 902, 929,
                                                                       2702, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5760, 0, 3, 6,
                                                                       5160, 2252, 5190, 5340,
                                                                       2396, 5370, 929, 956,
                                                                       2756, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5850, 0, 6, 5220,
                                                                       2324, 5250, 1010, 1028,
                                                                       2882, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5910, 0, 6, 5250,
                                                                       2342, 5280, 1028, 1046,
                                                                       2918, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5970, 0, 6, 5280,
                                                                       2360, 5310, 1046, 1064,
                                                                       2954, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6030, 0, 6, 5310,
                                                                       2378, 5340, 1064, 1082,
                                                                       2990, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6090, 0, 6, 5340,
                                                                       2396, 5370, 1082, 1100,
                                                                       3026, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 6150, 0, 3, 6,
                                                                       5400, 2540, 5490, 5850,
                                                                       2882, 5910, 1136, 1190,
                                                                       3278, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 6330, 0, 3, 6,
                                                                       5490, 2594, 5580, 5910,
                                                                       2918, 5970, 1190, 1244,
                                                                       3386, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 6510, 0, 3, 6,
                                                                       5580, 2648, 5670, 5970,
                                                                       2954, 6030, 1244, 1298,
                                                                       3494, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 6690, 0, 3, 6,
                                                                       5670, 2702, 5760, 6030,
                                                                       2990, 6090, 1298, 1352,
                                                                       3602, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6870, 0, 6, 5850,
                                                                       2882, 5910, 1460, 1490,
                                                                       3830, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6970, 0, 6, 5910,
                                                                       2918, 5970, 1490, 1520,
                                                                       3890, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7070, 0, 6, 5970,
                                                                       2954, 6030, 1520, 1550,
                                                                       3950, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7170, 0, 6, 6030,
                                                                       2990, 6090, 1550, 1580,
                                                                       4010, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 7270, 0, 3, 6,
                                                                       5400, 5490, 6150, 3278,
                                                                       6330, 6870, 3830, 6970,
                                                                       1640, 1730, 4430, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 7570, 0, 3, 6,
                                                                       5490, 5580, 6330, 3386,
                                                                       6510, 6970, 3890, 7070,
                                                                       1730, 1820, 4610, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 7870, 0, 3, 6,
                                                                       5580, 5670, 6510, 3494,
                                                                       6690, 7070, 3950, 7170,
                                                                       1820, 1910, 4790, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8170, 6, 2090,
                                                                       2096, 4970, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8185, 6, 2096,
                                                                       2102, 4980, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8200, 6, 2102,
                                                                       2108, 4990, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8215, 6, 2108,
                                                                       2114, 5000, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8230, 6, 2114,
                                                                       2120, 5010, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8245, 6, 2120,
                                                                       2126, 5020, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 8260, 6, 2126,
                                                                       2132, 5030, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8275, 3, 6, 8170,
                                                                       4970, 8185, 2144, 2162,
                                                                       5040, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8320, 3, 6, 8185,
                                                                       4980, 8200, 2162, 2180,
                                                                       5070, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8365, 3, 6, 8200,
                                                                       4990, 8215, 2180, 2198,
                                                                       5100, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8410, 3, 6, 8215,
                                                                       5000, 8230, 2198, 2216,
                                                                       5130, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8455, 3, 6, 8230,
                                                                       5010, 8245, 2216, 2234,
                                                                       5160, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8500, 3, 6, 8245,
                                                                       5020, 8260, 2234, 2252,
                                                                       5190, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8545, 0, 6, 8170,
                                                                       4970, 8185, 2288, 2306,
                                                                       5220, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8590, 0, 6, 8185,
                                                                       4980, 8200, 2306, 2324,
                                                                       5250, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8635, 0, 6, 8200,
                                                                       4990, 8215, 2324, 2342,
                                                                       5280, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8680, 0, 6, 8215,
                                                                       5000, 8230, 2342, 2360,
                                                                       5310, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8725, 0, 6, 8230,
                                                                       5010, 8245, 2360, 2378,
                                                                       5340, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8770, 0, 6, 8245,
                                                                       5020, 8260, 2378, 2396,
                                                                       5370, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 8815, 0, 3, 6,
                                                                       8275, 5040, 8320, 8545,
                                                                       5220, 8590, 2432, 2486,
                                                                       5400, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 8950, 0, 3, 6,
                                                                       8320, 5070, 8365, 8590,
                                                                       5250, 8635, 2486, 2540,
                                                                       5490, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9085, 0, 3, 6,
                                                                       8365, 5100, 8410, 8635,
                                                                       5280, 8680, 2540, 2594,
                                                                       5580, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9220, 0, 3, 6,
                                                                       8410, 5130, 8455, 8680,
                                                                       5310, 8725, 2594, 2648,
                                                                       5670, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9355, 0, 3, 6,
                                                                       8455, 5160, 8500, 8725,
                                                                       5340, 8770, 2648, 2702,
                                                                       5760, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 9490, 0, 6, 8545,
                                                                       5220, 8590, 2810, 2846,
                                                                       5850, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 9580, 0, 6, 8590,
                                                                       5250, 8635, 2846, 2882,
                                                                       5910, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 9670, 0, 6, 8635,
                                                                       5280, 8680, 2882, 2918,
                                                                       5970, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 9760, 0, 6, 8680,
                                                                       5310, 8725, 2918, 2954,
                                                                       6030, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 9850, 0, 6, 8725,
                                                                       5340, 8770, 2954, 2990,
                                                                       6090, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 9940, 0, 3, 6,
                                                                       8815, 5400, 8950, 9490,
                                                                       5850, 9580, 3062, 3170,
                                                                       6150, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 10210, 0, 3, 6,
                                                                       8950, 5490, 9085, 9580,
                                                                       5910, 9670, 3170, 3278,
                                                                       6330, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 10480, 0, 3, 6,
                                                                       9085, 5580, 9220, 9670,
                                                                       5970, 9760, 3278, 3386,
                                                                       6510, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 10750, 0, 3, 6,
                                                                       9220, 5670, 9355, 9760,
                                                                       6030, 9850, 3386, 3494,
                                                                       6690, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11020, 0, 6, 9490,
                                                                       5850, 9580, 3710, 3770,
                                                                       6870, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11170, 0, 6, 9580,
                                                                       5910, 9670, 3770, 3830,
                                                                       6970, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11320, 0, 6, 9670,
                                                                       5970, 9760, 3830, 3890,
                                                                       7070, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11470, 0, 6, 9760,
                                                                       6030, 9850, 3890, 3950,
                                                                       7170, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 11620, 0, 3, 6,
                                                                       8815, 8950, 9940, 6150,
                                                                       10210, 11020, 6870, 11170,
                                                                       4070, 4250, 7270, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 12070, 0, 3, 6,
                                                                       8950, 9085, 10210, 6330,
                                                                       10480, 11170, 6970, 11320,
                                                                       4250, 4430, 7570, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 12520, 0, 3, 6,
                                                                       9085, 9220, 10480, 6510,
                                                                       10750, 11320, 7070, 11470,
                                                                       4430, 4610, 7870, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12970, 6, 4970,
                                                                       4980, 8200, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12991, 6, 4980,
                                                                       4990, 8215, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13012, 6, 4990,
                                                                       5000, 8230, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13033, 6, 5000,
                                                                       5010, 8245, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 13054, 6, 5010,
                                                                       5020, 8260, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13075, 3, 6,
                                                                       12970, 8200, 12991, 5040,
                                                                       5070, 8365, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13138, 3, 6,
                                                                       12991, 8215, 13012, 5070,
                                                                       5100, 8410, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13201, 3, 6,
                                                                       13012, 8230, 13033, 5100,
                                                                       5130, 8455, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 13264, 3, 6,
                                                                       13033, 8245, 13054, 5130,
                                                                       5160, 8500, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13327, 0, 6,
                                                                       12970, 8200, 12991, 5220,
                                                                       5250, 8635, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13390, 0, 6,
                                                                       12991, 8215, 13012, 5250,
                                                                       5280, 8680, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13453, 0, 6,
                                                                       13012, 8230, 13033, 5280,
                                                                       5310, 8725, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13516, 0, 6,
                                                                       13033, 8245, 13054, 5310,
                                                                       5340, 8770, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 13579, 0, 3, 6,
                                                                       13075, 8365, 13138, 13327,
                                                                       8635, 13390, 5400, 5490,
                                                                       9085, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 13768, 0, 3, 6,
                                                                       13138, 8410, 13201, 13390,
                                                                       8680, 13453, 5490, 5580,
                                                                       9220, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 13957, 0, 3, 6,
                                                                       13201, 8455, 13264, 13453,
                                                                       8725, 13516, 5580, 5670,
                                                                       9355, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 14146, 0, 6,
                                                                       13327, 8635, 13390, 5850,
                                                                       5910, 9670, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 14272, 0, 6,
                                                                       13390, 8680, 13453, 5910,
                                                                       5970, 9760, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 14398, 0, 6,
                                                                       13453, 8725, 13516, 5970,
                                                                       6030, 9850, ncols, gamma,
                                                                       p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 14524, 0, 3, 6,
                                                                       13579, 9085, 13768, 14146,
                                                                       9670, 14272, 6150, 6330,
                                                                       10480, ncols, gamma, p,
                                                                       q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 14902, 0, 3, 6,
                                                                       13768, 9220, 13957, 14272,
                                                                       9760, 14398, 6330, 6510,
                                                                       10750, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 15280, 0, 6,
                                                                       14146, 9670, 14272, 6870,
                                                                       6970, 11320, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 15490, 0, 6,
                                                                       14272, 9760, 14398, 6970,
                                                                       7070, 11470, ncols, gamma,
                                                                       p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 15700, 0, 3, 6,
                                                                       13579, 13768, 14524,
                                                                       10480, 14902, 15280,
                                                                       11320, 15490, 7270, 7570,
                                                                       12520, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16330, 6, 8170,
                                                                       8185, 12970, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16358, 6, 8185,
                                                                       8200, 12991, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16386, 6, 8200,
                                                                       8215, 13012, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16414, 6, 8215,
                                                                       8230, 13033, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 16442, 6, 8230,
                                                                       8245, 13054, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 16470, 3, 6,
                                                                       16330, 12970, 16358, 8275,
                                                                       8320, 13075, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 16554, 3, 6,
                                                                       16358, 12991, 16386, 8320,
                                                                       8365, 13138, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 16638, 3, 6,
                                                                       16386, 13012, 16414, 8365,
                                                                       8410, 13201, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 16722, 3, 6,
                                                                       16414, 13033, 16442, 8410,
                                                                       8455, 13264, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16806, 0, 6,
                                                                       16330, 12970, 16358, 8545,
                                                                       8590, 13327, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16890, 0, 6,
                                                                       16358, 12991, 16386, 8590,
                                                                       8635, 13390, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16974, 0, 6,
                                                                       16386, 13012, 16414, 8635,
                                                                       8680, 13453, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 17058, 0, 6,
                                                                       16414, 13033, 16442, 8680,
                                                                       8725, 13516, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 17142, 0, 3, 6,
                                                                       16470, 13075, 16554,
                                                                       16806, 13327, 16890, 8815,
                                                                       8950, 13579, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 17394, 0, 3, 6,
                                                                       16554, 13138, 16638,
                                                                       16890, 13390, 16974, 8950,
                                                                       9085, 13768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 17646, 0, 3, 6,
                                                                       16638, 13201, 16722,
                                                                       16974, 13453, 17058, 9085,
                                                                       9220, 13957, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 17898, 0, 6,
                                                                       16806, 13327, 16890, 9490,
                                                                       9580, 14146, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 18066, 0, 6,
                                                                       16890, 13390, 16974, 9580,
                                                                       9670, 14272, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 18234, 0, 6,
                                                                       16974, 13453, 17058, 9670,
                                                                       9760, 14398, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 18402, 0, 3, 6,
                                                                       17142, 13579, 17394,
                                                                       17898, 14146, 18066, 9940,
                                                                       10210, 14524, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 18906, 0, 3, 6,
                                                                       17394, 13768, 17646,
                                                                       18066, 14272, 18234,
                                                                       10210, 10480, 14902,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 19410, 0, 6,
                                                                       17898, 14146, 18066,
                                                                       11020, 11170, 15280,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 19690, 0, 6,
                                                                       18066, 14272, 18234,
                                                                       11170, 11320, 15490,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpi_three_center_electron_repulsion_0(buffer, 19970, 0, 3, 6,
                                                                       17142, 17394, 18402,
                                                                       14524, 18906, 19410,
                                                                       15280, 19690, 11620,
                                                                       12070, 15700, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 20810, 19970, 10, 28, ncols, beta);

                    simdgeo::geom_s_y(buffer, 21090, 19970, 10, 28, ncols, beta);

                    simdgeo::geom_s_z(buffer, 21370, 19970, 10, 28, ncols, beta);

                    simdfunc::contract_primitives(buffer, 21650, 20810, 840, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 22490, 21650, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 22490, 13, nmax);

        simdtrf::transform_i_inner(buffer, 22490, 21930, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 91 * nvalues + n * npairs, nvalues, buffer, 22490,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 22490, 22210, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 182 * nvalues + n * npairs, nvalues, buffer, 22490,
                                   13, nmax);
    }

    for (size_t m = 0; m < 273; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
