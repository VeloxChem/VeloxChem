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


#include "SimdThreeCenterElectronRepulsionGeom010RecGSG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_gsg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_gsg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 15600, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 243 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 15600, 14790, 675, dimensions);

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

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 320, 0, 6, 47, 50,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 330, 0, 6, 50, 53,
                                                                       152, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 340, 0, 6, 53, 56,
                                                                       158, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 350, 0, 6, 56, 59,
                                                                       164, 170, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 360, 0, 6, 59, 62,
                                                                       170, 176, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 370, 0, 6, 62, 65,
                                                                       176, 182, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 380, 0, 6, 65, 68,
                                                                       182, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 390, 0, 3, 6, 146,
                                                                       152, 194, 212, 320, 330,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 420, 0, 3, 6, 152,
                                                                       158, 212, 230, 330, 340,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 450, 0, 3, 6, 158,
                                                                       164, 230, 248, 340, 350,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 480, 0, 3, 6, 164,
                                                                       170, 248, 266, 350, 360,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 510, 0, 3, 6, 170,
                                                                       176, 266, 284, 360, 370,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 540, 0, 3, 6, 176,
                                                                       182, 284, 302, 370, 380,
                                                                       ncols, gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 570, 0, 6, 146,
                                                                       152, 320, 330, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 585, 0, 6, 152,
                                                                       158, 330, 340, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 600, 0, 6, 158,
                                                                       164, 340, 350, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 615, 0, 6, 164,
                                                                       170, 350, 360, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 630, 0, 6, 170,
                                                                       176, 360, 370, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 645, 0, 6, 176,
                                                                       182, 370, 380, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 660, 0, 3, 6, 194,
                                                                       212, 320, 330, 390, 420,
                                                                       570, 585, ncols, gamma, p,
                                                                       q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 705, 0, 3, 6, 212,
                                                                       230, 330, 340, 420, 450,
                                                                       585, 600, ncols, gamma, p,
                                                                       q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 750, 0, 3, 6, 230,
                                                                       248, 340, 350, 450, 480,
                                                                       600, 615, ncols, gamma, p,
                                                                       q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 795, 0, 3, 6, 248,
                                                                       266, 350, 360, 480, 510,
                                                                       615, 630, ncols, gamma, p,
                                                                       q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 840, 0, 3, 6, 266,
                                                                       284, 360, 370, 510, 540,
                                                                       630, 645, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 885, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 888, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 891, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 894, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 897, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 900, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 903, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 906, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 909, 6, 12, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 918, 6, 13, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 927, 6, 14, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 936, 6, 15, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 945, 6, 16, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 954, 6, 17, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 963, 6, 18, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 972, 6, 12, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 981, 6, 13, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 990, 6, 14, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 999, 6, 15, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1008, 6, 16, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1017, 6, 17, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1026, 6, 18, 71,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1035, 6, 26, 53,
                                                                       92, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1062, 6, 29, 56,
                                                                       101, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1089, 6, 32, 59,
                                                                       110, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1116, 6, 35, 62,
                                                                       119, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1143, 6, 38, 65,
                                                                       128, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1170, 6, 41, 68,
                                                                       137, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1197, 6, 53, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1215, 6, 56, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1233, 6, 59, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1251, 6, 62, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1269, 6, 65, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1287, 6, 68, 188,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1305, 0, 6, 1035,
                                                                       92, 1062, 158, 230, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1359, 0, 6, 1062,
                                                                       101, 1089, 164, 248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1413, 0, 6, 1089,
                                                                       110, 1116, 170, 266,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1467, 0, 6, 1116,
                                                                       119, 1143, 176, 284,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1521, 0, 6, 1143,
                                                                       128, 1170, 182, 302,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1575, 6, 158, 340,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1605, 6, 164, 350,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1635, 6, 170, 360,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1665, 6, 176, 370,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1695, 6, 182, 380,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1725, 0, 6, 1305,
                                                                       230, 1359, 340, 450,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1815, 0, 6, 1359,
                                                                       248, 1413, 350, 480,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1905, 0, 6, 1413,
                                                                       266, 1467, 360, 510,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1995, 0, 6, 1467,
                                                                       284, 1521, 370, 540,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2085, 6, 340, 600,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2130, 6, 350, 615,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2175, 6, 360, 630,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2220, 6, 370, 645,
                                                                       ncols, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 2265, 0, 6, 1725,
                                                                       450, 1815, 600, 750,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 2400, 0, 6, 1815,
                                                                       480, 1905, 615, 795,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 2535, 0, 6, 1905,
                                                                       510, 1995, 630, 840,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2670, 6, 10, 11,
                                                                       885, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2676, 6, 11, 12,
                                                                       888, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2682, 6, 12, 13,
                                                                       891, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2688, 6, 13, 14,
                                                                       894, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2694, 6, 14, 15,
                                                                       897, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2700, 6, 15, 16,
                                                                       900, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2706, 6, 16, 17,
                                                                       903, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2712, 6, 17, 18,
                                                                       906, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2718, 3, 6, 2670,
                                                                       885, 2676, 909, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2736, 3, 6, 2676,
                                                                       888, 2682, 918, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2754, 3, 6, 2682,
                                                                       891, 2688, 927, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2772, 3, 6, 2688,
                                                                       894, 2694, 936, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2790, 3, 6, 2694,
                                                                       897, 2700, 945, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2808, 3, 6, 2700,
                                                                       900, 2706, 954, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2826, 3, 6, 2706,
                                                                       903, 2712, 963, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2844, 0, 6, 2670,
                                                                       885, 2676, 972, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2862, 0, 6, 2676,
                                                                       888, 2682, 981, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2880, 0, 6, 2682,
                                                                       891, 2688, 990, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2898, 0, 6, 2688,
                                                                       894, 2694, 999, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2916, 0, 6, 2694,
                                                                       897, 2700, 1008, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2934, 0, 6, 2700,
                                                                       900, 2706, 1017, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2952, 0, 6, 2706,
                                                                       903, 2712, 1026, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2970, 0, 3, 6,
                                                                       2718, 909, 2736, 2844,
                                                                       972, 2862, 74, 83, 1035,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3024, 0, 3, 6,
                                                                       2736, 918, 2754, 2862,
                                                                       981, 2880, 83, 92, 1062,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3078, 0, 3, 6,
                                                                       2754, 927, 2772, 2880,
                                                                       990, 2898, 92, 101, 1089,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3132, 0, 3, 6,
                                                                       2772, 936, 2790, 2898,
                                                                       999, 2916, 101, 110, 1116,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3186, 0, 3, 6,
                                                                       2790, 945, 2808, 2916,
                                                                       1008, 2934, 110, 119,
                                                                       1143, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3240, 0, 3, 6,
                                                                       2808, 954, 2826, 2934,
                                                                       1017, 2952, 119, 128,
                                                                       1170, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3294, 0, 6, 2844,
                                                                       972, 2862, 146, 152, 1197,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3330, 0, 6, 2862,
                                                                       981, 2880, 152, 158, 1215,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3366, 0, 6, 2880,
                                                                       990, 2898, 158, 164, 1233,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3402, 0, 6, 2898,
                                                                       999, 2916, 164, 170, 1251,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3438, 0, 6, 2916,
                                                                       1008, 2934, 170, 176,
                                                                       1269, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3474, 0, 6, 2934,
                                                                       1017, 2952, 176, 182,
                                                                       1287, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3510, 0, 3, 6,
                                                                       2970, 1035, 3024, 3294,
                                                                       1197, 3330, 194, 212,
                                                                       1305, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3618, 0, 3, 6,
                                                                       3024, 1062, 3078, 3330,
                                                                       1215, 3366, 212, 230,
                                                                       1359, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3726, 0, 3, 6,
                                                                       3078, 1089, 3132, 3366,
                                                                       1233, 3402, 230, 248,
                                                                       1413, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3834, 0, 3, 6,
                                                                       3132, 1116, 3186, 3402,
                                                                       1251, 3438, 248, 266,
                                                                       1467, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3942, 0, 3, 6,
                                                                       3186, 1143, 3240, 3438,
                                                                       1269, 3474, 266, 284,
                                                                       1521, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4050, 0, 6, 3294,
                                                                       1197, 3330, 320, 330,
                                                                       1575, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4110, 0, 6, 3330,
                                                                       1215, 3366, 330, 340,
                                                                       1605, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4170, 0, 6, 3366,
                                                                       1233, 3402, 340, 350,
                                                                       1635, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4230, 0, 6, 3402,
                                                                       1251, 3438, 350, 360,
                                                                       1665, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4290, 0, 6, 3438,
                                                                       1269, 3474, 360, 370,
                                                                       1695, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 4350, 0, 3, 6,
                                                                       2970, 3024, 3510, 1305,
                                                                       3618, 4050, 1575, 4110,
                                                                       390, 420, 1725, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 4530, 0, 3, 6,
                                                                       3024, 3078, 3618, 1359,
                                                                       3726, 4110, 1605, 4170,
                                                                       420, 450, 1815, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 4710, 0, 3, 6,
                                                                       3078, 3132, 3726, 1413,
                                                                       3834, 4170, 1635, 4230,
                                                                       450, 480, 1905, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 4890, 0, 3, 6,
                                                                       3132, 3186, 3834, 1467,
                                                                       3942, 4230, 1665, 4290,
                                                                       480, 510, 1995, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5070, 0, 6, 4050,
                                                                       1575, 4110, 570, 585,
                                                                       2085, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5160, 0, 6, 4110,
                                                                       1605, 4170, 585, 600,
                                                                       2130, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5250, 0, 6, 4170,
                                                                       1635, 4230, 600, 615,
                                                                       2175, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5340, 0, 6, 4230,
                                                                       1665, 4290, 615, 630,
                                                                       2220, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 5430, 0, 3, 6,
                                                                       3510, 3618, 4350, 1725,
                                                                       4530, 5070, 2085, 5160,
                                                                       660, 705, 2265, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 5700, 0, 3, 6,
                                                                       3618, 3726, 4530, 1815,
                                                                       4710, 5160, 2130, 5250,
                                                                       705, 750, 2400, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 5970, 0, 3, 6,
                                                                       3726, 3834, 4710, 1905,
                                                                       4890, 5250, 2175, 5340,
                                                                       750, 795, 2535, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6240, 6, 885, 888,
                                                                       2682, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6250, 6, 888, 891,
                                                                       2688, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6260, 6, 891, 894,
                                                                       2694, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6270, 6, 894, 897,
                                                                       2700, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6280, 6, 897, 900,
                                                                       2706, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6290, 6, 900, 903,
                                                                       2712, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6300, 3, 6, 6240,
                                                                       2682, 6250, 2754, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6330, 3, 6, 6250,
                                                                       2688, 6260, 2772, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6360, 3, 6, 6260,
                                                                       2694, 6270, 2790, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6390, 3, 6, 6270,
                                                                       2700, 6280, 2808, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6420, 3, 6, 6280,
                                                                       2706, 6290, 2826, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6450, 0, 6, 6240,
                                                                       2682, 6250, 972, 981,
                                                                       2880, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6480, 0, 6, 6250,
                                                                       2688, 6260, 981, 990,
                                                                       2898, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6510, 0, 6, 6260,
                                                                       2694, 6270, 990, 999,
                                                                       2916, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6540, 0, 6, 6270,
                                                                       2700, 6280, 999, 1008,
                                                                       2934, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 6570, 0, 6, 6280,
                                                                       2706, 6290, 1008, 1017,
                                                                       2952, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6600, 0, 3, 6,
                                                                       6300, 2754, 6330, 6450,
                                                                       2880, 6480, 1035, 1062,
                                                                       3078, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6690, 0, 3, 6,
                                                                       6330, 2772, 6360, 6480,
                                                                       2898, 6510, 1062, 1089,
                                                                       3132, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6780, 0, 3, 6,
                                                                       6360, 2790, 6390, 6510,
                                                                       2916, 6540, 1089, 1116,
                                                                       3186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6870, 0, 3, 6,
                                                                       6390, 2808, 6420, 6540,
                                                                       2934, 6570, 1116, 1143,
                                                                       3240, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6960, 0, 6, 6450,
                                                                       2880, 6480, 1197, 1215,
                                                                       3366, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7020, 0, 6, 6480,
                                                                       2898, 6510, 1215, 1233,
                                                                       3402, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7080, 0, 6, 6510,
                                                                       2916, 6540, 1233, 1251,
                                                                       3438, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7140, 0, 6, 6540,
                                                                       2934, 6570, 1251, 1269,
                                                                       3474, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 7200, 0, 3, 6,
                                                                       6600, 3078, 6690, 6960,
                                                                       3366, 7020, 1305, 1359,
                                                                       3726, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 7380, 0, 3, 6,
                                                                       6690, 3132, 6780, 7020,
                                                                       3402, 7080, 1359, 1413,
                                                                       3834, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 7560, 0, 3, 6,
                                                                       6780, 3186, 6870, 7080,
                                                                       3438, 7140, 1413, 1467,
                                                                       3942, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7740, 0, 6, 6960,
                                                                       3366, 7020, 1575, 1605,
                                                                       4170, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7840, 0, 6, 7020,
                                                                       3402, 7080, 1605, 1635,
                                                                       4230, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 7940, 0, 6, 7080,
                                                                       3438, 7140, 1635, 1665,
                                                                       4290, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 8040, 0, 3, 6,
                                                                       6600, 6690, 7200, 3726,
                                                                       7380, 7740, 4170, 7840,
                                                                       1725, 1815, 4710, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 8340, 0, 3, 6,
                                                                       6690, 6780, 7380, 3834,
                                                                       7560, 7840, 4230, 7940,
                                                                       1815, 1905, 4890, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8640, 0, 6, 7740,
                                                                       4170, 7840, 2085, 2130,
                                                                       5250, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 8790, 0, 6, 7840,
                                                                       4230, 7940, 2130, 2175,
                                                                       5340, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 8940, 0, 3, 6,
                                                                       7200, 7380, 8040, 4710,
                                                                       8340, 8640, 5250, 8790,
                                                                       2265, 2400, 5970, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9390, 6, 2670,
                                                                       2676, 6240, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9405, 6, 2676,
                                                                       2682, 6250, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9420, 6, 2682,
                                                                       2688, 6260, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9435, 6, 2688,
                                                                       2694, 6270, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9450, 6, 2694,
                                                                       2700, 6280, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9465, 6, 2700,
                                                                       2706, 6290, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9480, 3, 6, 9390,
                                                                       6240, 9405, 2718, 2736,
                                                                       6300, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9525, 3, 6, 9405,
                                                                       6250, 9420, 2736, 2754,
                                                                       6330, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9570, 3, 6, 9420,
                                                                       6260, 9435, 2754, 2772,
                                                                       6360, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9615, 3, 6, 9435,
                                                                       6270, 9450, 2772, 2790,
                                                                       6390, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 9660, 3, 6, 9450,
                                                                       6280, 9465, 2790, 2808,
                                                                       6420, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9705, 0, 6, 9390,
                                                                       6240, 9405, 2844, 2862,
                                                                       6450, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9750, 0, 6, 9405,
                                                                       6250, 9420, 2862, 2880,
                                                                       6480, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9795, 0, 6, 9420,
                                                                       6260, 9435, 2880, 2898,
                                                                       6510, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9840, 0, 6, 9435,
                                                                       6270, 9450, 2898, 2916,
                                                                       6540, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9885, 0, 6, 9450,
                                                                       6280, 9465, 2916, 2934,
                                                                       6570, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9930, 0, 3, 6,
                                                                       9480, 6300, 9525, 9705,
                                                                       6450, 9750, 2970, 3024,
                                                                       6600, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 10065, 0, 3, 6,
                                                                       9525, 6330, 9570, 9750,
                                                                       6480, 9795, 3024, 3078,
                                                                       6690, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 10200, 0, 3, 6,
                                                                       9570, 6360, 9615, 9795,
                                                                       6510, 9840, 3078, 3132,
                                                                       6780, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 10335, 0, 3, 6,
                                                                       9615, 6390, 9660, 9840,
                                                                       6540, 9885, 3132, 3186,
                                                                       6870, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10470, 0, 6, 9705,
                                                                       6450, 9750, 3294, 3330,
                                                                       6960, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10560, 0, 6, 9750,
                                                                       6480, 9795, 3330, 3366,
                                                                       7020, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10650, 0, 6, 9795,
                                                                       6510, 9840, 3366, 3402,
                                                                       7080, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10740, 0, 6, 9840,
                                                                       6540, 9885, 3402, 3438,
                                                                       7140, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 10830, 0, 3, 6,
                                                                       9930, 6600, 10065, 10470,
                                                                       6960, 10560, 3510, 3618,
                                                                       7200, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 11100, 0, 3, 6,
                                                                       10065, 6690, 10200, 10560,
                                                                       7020, 10650, 3618, 3726,
                                                                       7380, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 11370, 0, 3, 6,
                                                                       10200, 6780, 10335, 10650,
                                                                       7080, 10740, 3726, 3834,
                                                                       7560, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11640, 0, 6,
                                                                       10470, 6960, 10560, 4050,
                                                                       4110, 7740, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11790, 0, 6,
                                                                       10560, 7020, 10650, 4110,
                                                                       4170, 7840, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 11940, 0, 6,
                                                                       10650, 7080, 10740, 4170,
                                                                       4230, 7940, ncols, gamma,
                                                                       p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 12090, 0, 3, 6,
                                                                       9930, 10065, 10830, 7200,
                                                                       11100, 11640, 7740, 11790,
                                                                       4350, 4530, 8040, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 12540, 0, 3, 6,
                                                                       10065, 10200, 11100, 7380,
                                                                       11370, 11790, 7840, 11940,
                                                                       4530, 4710, 8340, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 12990, 0, 6,
                                                                       11640, 7740, 11790, 5070,
                                                                       5160, 8640, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 13215, 0, 6,
                                                                       11790, 7840, 11940, 5160,
                                                                       5250, 8790, ncols, gamma,
                                                                       p, q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 13440, 0, 3, 6,
                                                                       10830, 11100, 12090, 8040,
                                                                       12540, 12990, 8640, 13215,
                                                                       5430, 5700, 8940, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 14115, 13440, 15, 15, ncols, beta);

                    simdgeo::geom_s_y(buffer, 14340, 13440, 15, 15, ncols, beta);

                    simdgeo::geom_s_z(buffer, 14565, 13440, 15, 15, ncols, beta);

                    simdfunc::contract_primitives(buffer, 14790, 14115, 675, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 15465, 14790, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 15465, 9, nmax);

        simdtrf::transform_g_inner(buffer, 15465, 15015, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 81 * nvalues + n * npairs, nvalues, buffer, 15465, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 15465, 15240, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 162 * nvalues + n * npairs, nvalues, buffer, 15465,
                                   9, nmax);
    }

    for (size_t m = 0; m < 243; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
