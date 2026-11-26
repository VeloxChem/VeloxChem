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


#include "SimdThreeCenterElectronRepulsionGeom010RecGSH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_gsh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_gsh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 26407, 0, 0, dimensions);

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

    simdfunc::compute_pair_exponents(a_function, b_function, coordinates, nmax);

    for (size_t n = 0; n < natoms; n++)
    {
        simdfunc::prepare_buffer(buffer, 26407, 25297, 945, dimensions);

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
                                                        5, 6, 7, 8, 9, 10}, ncols, fj,
                                                        i * nprim_b + j, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 885, 6, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 888, 6, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 891, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 894, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 897, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 900, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 903, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 906, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 909, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 912, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 915, 6, 12, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 924, 6, 13, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 933, 6, 14, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 942, 6, 15, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 951, 6, 16, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 960, 6, 17, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 969, 6, 18, 44,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 978, 6, 10, 47,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 987, 6, 11, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 996, 6, 12, 53,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1005, 6, 13, 56,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1014, 6, 14, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1023, 6, 15, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1032, 6, 16, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1041, 6, 17, 68,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1050, 6, 18, 71,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1059, 6, 20, 47,
                                                                       74, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1086, 6, 23, 50,
                                                                       83, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1113, 6, 26, 53,
                                                                       92, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1140, 6, 29, 56,
                                                                       101, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1167, 6, 32, 59,
                                                                       110, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1194, 6, 35, 62,
                                                                       119, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1221, 6, 38, 65,
                                                                       128, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1248, 6, 41, 68,
                                                                       137, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1275, 6, 47, 146,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1293, 6, 50, 152,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1311, 6, 53, 158,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1329, 6, 56, 164,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1347, 6, 59, 170,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1365, 6, 62, 176,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1383, 6, 65, 182,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1401, 6, 68, 188,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1419, 0, 6, 1059,
                                                                       74, 1086, 146, 194, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1473, 0, 6, 1086,
                                                                       83, 1113, 152, 212, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1527, 0, 6, 1113,
                                                                       92, 1140, 158, 230, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1581, 0, 6, 1140,
                                                                       101, 1167, 164, 248,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1635, 0, 6, 1167,
                                                                       110, 1194, 170, 266,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1689, 0, 6, 1194,
                                                                       119, 1221, 176, 284,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1743, 0, 6, 1221,
                                                                       128, 1248, 182, 302,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1797, 6, 146, 320,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1827, 6, 152, 330,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1857, 6, 158, 340,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1887, 6, 164, 350,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1917, 6, 170, 360,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1947, 6, 176, 370,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1977, 6, 182, 380,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2007, 0, 6, 1419,
                                                                       194, 1473, 320, 390,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2097, 0, 6, 1473,
                                                                       212, 1527, 330, 420,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2187, 0, 6, 1527,
                                                                       230, 1581, 340, 450,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2277, 0, 6, 1581,
                                                                       248, 1635, 350, 480,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2367, 0, 6, 1635,
                                                                       266, 1689, 360, 510,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2457, 0, 6, 1689,
                                                                       284, 1743, 370, 540,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2547, 6, 320, 570,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2592, 6, 330, 585,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2637, 6, 340, 600,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2682, 6, 350, 615,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2727, 6, 360, 630,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2772, 6, 370, 645,
                                                                       ncols, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 2817, 0, 6, 2007,
                                                                       390, 2097, 570, 660,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 2952, 0, 6, 2097,
                                                                       420, 2187, 585, 705,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 3087, 0, 6, 2187,
                                                                       450, 2277, 600, 750,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 3222, 0, 6, 2277,
                                                                       480, 2367, 615, 795,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 3357, 0, 6, 2367,
                                                                       510, 2457, 630, 840,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3492, 6, 10, 11,
                                                                       891, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3498, 6, 11, 12,
                                                                       894, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3504, 6, 12, 13,
                                                                       897, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3510, 6, 13, 14,
                                                                       900, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3516, 6, 14, 15,
                                                                       903, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3522, 6, 15, 16,
                                                                       906, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3528, 6, 16, 17,
                                                                       909, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3534, 6, 17, 18,
                                                                       912, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3540, 3, 6, 3492,
                                                                       891, 3498, 915, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3558, 3, 6, 3498,
                                                                       894, 3504, 924, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3576, 3, 6, 3504,
                                                                       897, 3510, 933, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3594, 3, 6, 3510,
                                                                       900, 3516, 942, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3612, 3, 6, 3516,
                                                                       903, 3522, 951, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3630, 3, 6, 3522,
                                                                       906, 3528, 960, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3648, 3, 6, 3528,
                                                                       909, 3534, 969, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3666, 0, 6, 3492,
                                                                       891, 3498, 996, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3684, 0, 6, 3498,
                                                                       894, 3504, 1005, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3702, 0, 6, 3504,
                                                                       897, 3510, 1014, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3720, 0, 6, 3510,
                                                                       900, 3516, 1023, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3738, 0, 6, 3516,
                                                                       903, 3522, 1032, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3756, 0, 6, 3522,
                                                                       906, 3528, 1041, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3774, 0, 6, 3528,
                                                                       909, 3534, 1050, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3792, 0, 3, 6,
                                                                       3540, 915, 3558, 3666,
                                                                       996, 3684, 74, 83, 1113,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3846, 0, 3, 6,
                                                                       3558, 924, 3576, 3684,
                                                                       1005, 3702, 83, 92, 1140,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3900, 0, 3, 6,
                                                                       3576, 933, 3594, 3702,
                                                                       1014, 3720, 92, 101, 1167,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3954, 0, 3, 6,
                                                                       3594, 942, 3612, 3720,
                                                                       1023, 3738, 101, 110,
                                                                       1194, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4008, 0, 3, 6,
                                                                       3612, 951, 3630, 3738,
                                                                       1032, 3756, 110, 119,
                                                                       1221, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4062, 0, 3, 6,
                                                                       3630, 960, 3648, 3756,
                                                                       1041, 3774, 119, 128,
                                                                       1248, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4116, 0, 6, 3666,
                                                                       996, 3684, 146, 152, 1311,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4152, 0, 6, 3684,
                                                                       1005, 3702, 152, 158,
                                                                       1329, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4188, 0, 6, 3702,
                                                                       1014, 3720, 158, 164,
                                                                       1347, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4224, 0, 6, 3720,
                                                                       1023, 3738, 164, 170,
                                                                       1365, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4260, 0, 6, 3738,
                                                                       1032, 3756, 170, 176,
                                                                       1383, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4296, 0, 6, 3756,
                                                                       1041, 3774, 176, 182,
                                                                       1401, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4332, 0, 3, 6,
                                                                       3792, 1113, 3846, 4116,
                                                                       1311, 4152, 194, 212,
                                                                       1527, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4440, 0, 3, 6,
                                                                       3846, 1140, 3900, 4152,
                                                                       1329, 4188, 212, 230,
                                                                       1581, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4548, 0, 3, 6,
                                                                       3900, 1167, 3954, 4188,
                                                                       1347, 4224, 230, 248,
                                                                       1635, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4656, 0, 3, 6,
                                                                       3954, 1194, 4008, 4224,
                                                                       1365, 4260, 248, 266,
                                                                       1689, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4764, 0, 3, 6,
                                                                       4008, 1221, 4062, 4260,
                                                                       1383, 4296, 266, 284,
                                                                       1743, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4872, 0, 6, 4116,
                                                                       1311, 4152, 320, 330,
                                                                       1857, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4932, 0, 6, 4152,
                                                                       1329, 4188, 330, 340,
                                                                       1887, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4992, 0, 6, 4188,
                                                                       1347, 4224, 340, 350,
                                                                       1917, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5052, 0, 6, 4224,
                                                                       1365, 4260, 350, 360,
                                                                       1947, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5112, 0, 6, 4260,
                                                                       1383, 4296, 360, 370,
                                                                       1977, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 5172, 0, 3, 6,
                                                                       3792, 3846, 4332, 1527,
                                                                       4440, 4872, 1857, 4932,
                                                                       390, 420, 2187, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 5352, 0, 3, 6,
                                                                       3846, 3900, 4440, 1581,
                                                                       4548, 4932, 1887, 4992,
                                                                       420, 450, 2277, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 5532, 0, 3, 6,
                                                                       3900, 3954, 4548, 1635,
                                                                       4656, 4992, 1917, 5052,
                                                                       450, 480, 2367, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 5712, 0, 3, 6,
                                                                       3954, 4008, 4656, 1689,
                                                                       4764, 5052, 1947, 5112,
                                                                       480, 510, 2457, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5892, 0, 6, 4872,
                                                                       1857, 4932, 570, 585,
                                                                       2637, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5982, 0, 6, 4932,
                                                                       1887, 4992, 585, 600,
                                                                       2682, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6072, 0, 6, 4992,
                                                                       1917, 5052, 600, 615,
                                                                       2727, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6162, 0, 6, 5052,
                                                                       1947, 5112, 615, 630,
                                                                       2772, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 6252, 0, 3, 6,
                                                                       4332, 4440, 5172, 2187,
                                                                       5352, 5892, 2637, 5982,
                                                                       660, 705, 3087, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 6522, 0, 3, 6,
                                                                       4440, 4548, 5352, 2277,
                                                                       5532, 5982, 2682, 6072,
                                                                       705, 750, 3222, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 6792, 0, 3, 6,
                                                                       4548, 4656, 5532, 2367,
                                                                       5712, 6072, 2727, 6162,
                                                                       750, 795, 3357, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7062, 6, 885, 888,
                                                                       3492, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7072, 6, 888, 891,
                                                                       3498, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7082, 6, 891, 894,
                                                                       3504, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7092, 6, 894, 897,
                                                                       3510, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7102, 6, 897, 900,
                                                                       3516, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7112, 6, 900, 903,
                                                                       3522, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7122, 6, 903, 906,
                                                                       3528, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7132, 6, 906, 909,
                                                                       3534, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7142, 3, 6, 7062,
                                                                       3492, 7072, 3540, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7172, 3, 6, 7072,
                                                                       3498, 7082, 3558, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7202, 3, 6, 7082,
                                                                       3504, 7092, 3576, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7232, 3, 6, 7092,
                                                                       3510, 7102, 3594, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7262, 3, 6, 7102,
                                                                       3516, 7112, 3612, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7292, 3, 6, 7112,
                                                                       3522, 7122, 3630, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7322, 3, 6, 7122,
                                                                       3528, 7132, 3648, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7352, 0, 6, 7062,
                                                                       3492, 7072, 978, 987,
                                                                       3666, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7382, 0, 6, 7072,
                                                                       3498, 7082, 987, 996,
                                                                       3684, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7412, 0, 6, 7082,
                                                                       3504, 7092, 996, 1005,
                                                                       3702, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7442, 0, 6, 7092,
                                                                       3510, 7102, 1005, 1014,
                                                                       3720, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7472, 0, 6, 7102,
                                                                       3516, 7112, 1014, 1023,
                                                                       3738, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7502, 0, 6, 7112,
                                                                       3522, 7122, 1023, 1032,
                                                                       3756, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7532, 0, 6, 7122,
                                                                       3528, 7132, 1032, 1041,
                                                                       3774, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7562, 0, 3, 6,
                                                                       7142, 3540, 7172, 7352,
                                                                       3666, 7382, 1059, 1086,
                                                                       3792, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7652, 0, 3, 6,
                                                                       7172, 3558, 7202, 7382,
                                                                       3684, 7412, 1086, 1113,
                                                                       3846, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7742, 0, 3, 6,
                                                                       7202, 3576, 7232, 7412,
                                                                       3702, 7442, 1113, 1140,
                                                                       3900, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7832, 0, 3, 6,
                                                                       7232, 3594, 7262, 7442,
                                                                       3720, 7472, 1140, 1167,
                                                                       3954, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7922, 0, 3, 6,
                                                                       7262, 3612, 7292, 7472,
                                                                       3738, 7502, 1167, 1194,
                                                                       4008, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 8012, 0, 3, 6,
                                                                       7292, 3630, 7322, 7502,
                                                                       3756, 7532, 1194, 1221,
                                                                       4062, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8102, 0, 6, 7352,
                                                                       3666, 7382, 1275, 1293,
                                                                       4116, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8162, 0, 6, 7382,
                                                                       3684, 7412, 1293, 1311,
                                                                       4152, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8222, 0, 6, 7412,
                                                                       3702, 7442, 1311, 1329,
                                                                       4188, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8282, 0, 6, 7442,
                                                                       3720, 7472, 1329, 1347,
                                                                       4224, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8342, 0, 6, 7472,
                                                                       3738, 7502, 1347, 1365,
                                                                       4260, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8402, 0, 6, 7502,
                                                                       3756, 7532, 1365, 1383,
                                                                       4296, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 8462, 0, 3, 6,
                                                                       7562, 3792, 7652, 8102,
                                                                       4116, 8162, 1419, 1473,
                                                                       4332, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 8642, 0, 3, 6,
                                                                       7652, 3846, 7742, 8162,
                                                                       4152, 8222, 1473, 1527,
                                                                       4440, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 8822, 0, 3, 6,
                                                                       7742, 3900, 7832, 8222,
                                                                       4188, 8282, 1527, 1581,
                                                                       4548, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 9002, 0, 3, 6,
                                                                       7832, 3954, 7922, 8282,
                                                                       4224, 8342, 1581, 1635,
                                                                       4656, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 9182, 0, 3, 6,
                                                                       7922, 4008, 8012, 8342,
                                                                       4260, 8402, 1635, 1689,
                                                                       4764, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9362, 0, 6, 8102,
                                                                       4116, 8162, 1797, 1827,
                                                                       4872, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9462, 0, 6, 8162,
                                                                       4152, 8222, 1827, 1857,
                                                                       4932, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9562, 0, 6, 8222,
                                                                       4188, 8282, 1857, 1887,
                                                                       4992, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9662, 0, 6, 8282,
                                                                       4224, 8342, 1887, 1917,
                                                                       5052, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9762, 0, 6, 8342,
                                                                       4260, 8402, 1917, 1947,
                                                                       5112, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 9862, 0, 3, 6,
                                                                       7562, 7652, 8462, 4332,
                                                                       8642, 9362, 4872, 9462,
                                                                       2007, 2097, 5172, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 10162, 0, 3, 6,
                                                                       7652, 7742, 8642, 4440,
                                                                       8822, 9462, 4932, 9562,
                                                                       2097, 2187, 5352, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 10462, 0, 3, 6,
                                                                       7742, 7832, 8822, 4548,
                                                                       9002, 9562, 4992, 9662,
                                                                       2187, 2277, 5532, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 10762, 0, 3, 6,
                                                                       7832, 7922, 9002, 4656,
                                                                       9182, 9662, 5052, 9762,
                                                                       2277, 2367, 5712, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11062, 0, 6, 9362,
                                                                       4872, 9462, 2547, 2592,
                                                                       5892, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11212, 0, 6, 9462,
                                                                       4932, 9562, 2592, 2637,
                                                                       5982, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11362, 0, 6, 9562,
                                                                       4992, 9662, 2637, 2682,
                                                                       6072, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11512, 0, 6, 9662,
                                                                       5052, 9762, 2682, 2727,
                                                                       6162, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 11662, 0, 3, 6,
                                                                       8462, 8642, 9862, 5172,
                                                                       10162, 11062, 5892, 11212,
                                                                       2817, 2952, 6252, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 12112, 0, 3, 6,
                                                                       8642, 8822, 10162, 5352,
                                                                       10462, 11212, 5982, 11362,
                                                                       2952, 3087, 6522, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 12562, 0, 3, 6,
                                                                       8822, 9002, 10462, 5532,
                                                                       10762, 11362, 6072, 11512,
                                                                       3087, 3222, 6792, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13012, 6, 3492,
                                                                       3498, 7082, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13027, 6, 3498,
                                                                       3504, 7092, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13042, 6, 3504,
                                                                       3510, 7102, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13057, 6, 3510,
                                                                       3516, 7112, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13072, 6, 3516,
                                                                       3522, 7122, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 13087, 6, 3522,
                                                                       3528, 7132, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13102, 3, 6,
                                                                       13012, 7082, 13027, 3540,
                                                                       3558, 7202, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13147, 3, 6,
                                                                       13027, 7092, 13042, 3558,
                                                                       3576, 7232, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13192, 3, 6,
                                                                       13042, 7102, 13057, 3576,
                                                                       3594, 7262, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13237, 3, 6,
                                                                       13057, 7112, 13072, 3594,
                                                                       3612, 7292, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 13282, 3, 6,
                                                                       13072, 7122, 13087, 3612,
                                                                       3630, 7322, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 13327, 0, 6,
                                                                       13012, 7082, 13027, 3666,
                                                                       3684, 7412, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 13372, 0, 6,
                                                                       13027, 7092, 13042, 3684,
                                                                       3702, 7442, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 13417, 0, 6,
                                                                       13042, 7102, 13057, 3702,
                                                                       3720, 7472, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 13462, 0, 6,
                                                                       13057, 7112, 13072, 3720,
                                                                       3738, 7502, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 13507, 0, 6,
                                                                       13072, 7122, 13087, 3738,
                                                                       3756, 7532, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 13552, 0, 3, 6,
                                                                       13102, 7202, 13147, 13327,
                                                                       7412, 13372, 3792, 3846,
                                                                       7742, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 13687, 0, 3, 6,
                                                                       13147, 7232, 13192, 13372,
                                                                       7442, 13417, 3846, 3900,
                                                                       7832, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 13822, 0, 3, 6,
                                                                       13192, 7262, 13237, 13417,
                                                                       7472, 13462, 3900, 3954,
                                                                       7922, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 13957, 0, 3, 6,
                                                                       13237, 7292, 13282, 13462,
                                                                       7502, 13507, 3954, 4008,
                                                                       8012, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14092, 0, 6,
                                                                       13327, 7412, 13372, 4116,
                                                                       4152, 8222, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14182, 0, 6,
                                                                       13372, 7442, 13417, 4152,
                                                                       4188, 8282, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14272, 0, 6,
                                                                       13417, 7472, 13462, 4188,
                                                                       4224, 8342, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 14362, 0, 6,
                                                                       13462, 7502, 13507, 4224,
                                                                       4260, 8402, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 14452, 0, 3, 6,
                                                                       13552, 7742, 13687, 14092,
                                                                       8222, 14182, 4332, 4440,
                                                                       8822, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 14722, 0, 3, 6,
                                                                       13687, 7832, 13822, 14182,
                                                                       8282, 14272, 4440, 4548,
                                                                       9002, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 14992, 0, 3, 6,
                                                                       13822, 7922, 13957, 14272,
                                                                       8342, 14362, 4548, 4656,
                                                                       9182, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15262, 0, 6,
                                                                       14092, 8222, 14182, 4872,
                                                                       4932, 9562, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15412, 0, 6,
                                                                       14182, 8282, 14272, 4932,
                                                                       4992, 9662, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 15562, 0, 6,
                                                                       14272, 8342, 14362, 4992,
                                                                       5052, 9762, ncols, gamma,
                                                                       p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 15712, 0, 3, 6,
                                                                       13552, 13687, 14452, 8822,
                                                                       14722, 15262, 9562, 15412,
                                                                       5172, 5352, 10462, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 16162, 0, 3, 6,
                                                                       13687, 13822, 14722, 9002,
                                                                       14992, 15412, 9662, 15562,
                                                                       5352, 5532, 10762, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 16612, 0, 6,
                                                                       15262, 9562, 15412, 5892,
                                                                       5982, 11362, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 16837, 0, 6,
                                                                       15412, 9662, 15562, 5982,
                                                                       6072, 11512, ncols, gamma,
                                                                       p, q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 17062, 0, 3, 6,
                                                                       14452, 14722, 15712,
                                                                       10462, 16162, 16612,
                                                                       11362, 16837, 6252, 6522,
                                                                       12562, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17737, 6, 7062,
                                                                       7072, 13012, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17758, 6, 7072,
                                                                       7082, 13027, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17779, 6, 7082,
                                                                       7092, 13042, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17800, 6, 7092,
                                                                       7102, 13057, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17821, 6, 7102,
                                                                       7112, 13072, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 17842, 6, 7112,
                                                                       7122, 13087, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17863, 3, 6,
                                                                       17737, 13012, 17758, 7142,
                                                                       7172, 13102, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17926, 3, 6,
                                                                       17758, 13027, 17779, 7172,
                                                                       7202, 13147, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 17989, 3, 6,
                                                                       17779, 13042, 17800, 7202,
                                                                       7232, 13192, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 18052, 3, 6,
                                                                       17800, 13057, 17821, 7232,
                                                                       7262, 13237, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 18115, 3, 6,
                                                                       17821, 13072, 17842, 7262,
                                                                       7292, 13282, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 18178, 0, 6,
                                                                       17737, 13012, 17758, 7352,
                                                                       7382, 13327, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 18241, 0, 6,
                                                                       17758, 13027, 17779, 7382,
                                                                       7412, 13372, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 18304, 0, 6,
                                                                       17779, 13042, 17800, 7412,
                                                                       7442, 13417, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 18367, 0, 6,
                                                                       17800, 13057, 17821, 7442,
                                                                       7472, 13462, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 18430, 0, 6,
                                                                       17821, 13072, 17842, 7472,
                                                                       7502, 13507, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 18493, 0, 3, 6,
                                                                       17863, 13102, 17926,
                                                                       18178, 13327, 18241, 7562,
                                                                       7652, 13552, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 18682, 0, 3, 6,
                                                                       17926, 13147, 17989,
                                                                       18241, 13372, 18304, 7652,
                                                                       7742, 13687, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 18871, 0, 3, 6,
                                                                       17989, 13192, 18052,
                                                                       18304, 13417, 18367, 7742,
                                                                       7832, 13822, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 19060, 0, 3, 6,
                                                                       18052, 13237, 18115,
                                                                       18367, 13462, 18430, 7832,
                                                                       7922, 13957, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 19249, 0, 6,
                                                                       18178, 13327, 18241, 8102,
                                                                       8162, 14092, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 19375, 0, 6,
                                                                       18241, 13372, 18304, 8162,
                                                                       8222, 14182, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 19501, 0, 6,
                                                                       18304, 13417, 18367, 8222,
                                                                       8282, 14272, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 19627, 0, 6,
                                                                       18367, 13462, 18430, 8282,
                                                                       8342, 14362, ncols, gamma,
                                                                       p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 19753, 0, 3, 6,
                                                                       18493, 13552, 18682,
                                                                       19249, 14092, 19375, 8462,
                                                                       8642, 14452, ncols, gamma,
                                                                       p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 20131, 0, 3, 6,
                                                                       18682, 13687, 18871,
                                                                       19375, 14182, 19501, 8642,
                                                                       8822, 14722, ncols, gamma,
                                                                       p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 20509, 0, 3, 6,
                                                                       18871, 13822, 19060,
                                                                       19501, 14272, 19627, 8822,
                                                                       9002, 14992, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 20887, 0, 6,
                                                                       19249, 14092, 19375, 9362,
                                                                       9462, 15262, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 21097, 0, 6,
                                                                       19375, 14182, 19501, 9462,
                                                                       9562, 15412, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 21307, 0, 6,
                                                                       19501, 14272, 19627, 9562,
                                                                       9662, 15562, ncols, gamma,
                                                                       p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 21517, 0, 3, 6,
                                                                       18493, 18682, 19753,
                                                                       14452, 20131, 20887,
                                                                       15262, 21097, 9862, 10162,
                                                                       15712, ncols, gamma, p,
                                                                       q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 22147, 0, 3, 6,
                                                                       18682, 18871, 20131,
                                                                       14722, 20509, 21097,
                                                                       15412, 21307, 10162,
                                                                       10462, 16162, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 22777, 0, 6,
                                                                       20887, 15262, 21097,
                                                                       11062, 11212, 16612,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 23092, 0, 6,
                                                                       21097, 15412, 21307,
                                                                       11212, 11362, 16837,
                                                                       ncols, gamma, p, q);

                    compute_prim_gph_three_center_electron_repulsion_0(buffer, 23407, 0, 3, 6,
                                                                       19753, 20131, 21517,
                                                                       15712, 22147, 22777,
                                                                       16612, 23092, 11662,
                                                                       12112, 17062, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 24352, 23407, 15, 21, ncols, beta);

                    simdgeo::geom_s_y(buffer, 24667, 23407, 15, 21, ncols, beta);

                    simdgeo::geom_s_z(buffer, 24982, 23407, 15, 21, ncols, beta);

                    simdfunc::contract_primitives(buffer, 25297, 24352, 945, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 26242, 25297, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 26242, 11, nmax);

        simdtrf::transform_h_inner(buffer, 26242, 25612, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 99 * nvalues + n * npairs, nvalues, buffer, 26242,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 26242, 25927, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 198 * nvalues + n * npairs, nvalues, buffer, 26242,
                                   11, nmax);
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
