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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecPSI.hpp"

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
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_psi_three_center_electron_repulsion(double               *values,
                                                        const size_t          npairs,
                                                        const size_t          natoms,
                                                        const CBasisFunction &a_function,
                                                        const CBasisFunction &b_function,
                                                        const CBasisFunction &c_function,
                                                        const CSimdMatrix    &coordinates,
                                                        const CSimdMatrix    &c_coordinates,
                                                        CSimdMatrix          &buffer,
                                                        const double          omega,
                                                        const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_geom_010_psi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 8022, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 234 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 8022, 7479, 504, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 9, 6, 8,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 19, 6, 8,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 3, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 56, 3, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 59, 3, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 62, 3, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 65, 3, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 68, 3, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 71, 3, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 74, 3, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 77, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 80, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 83, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 86, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 89, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 92, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 95, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 98, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 101, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 104, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 107, 0, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 110, 0, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 113, 0, 6, 10, 11,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 122, 0, 6, 11, 12,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 131, 0, 6, 12, 13,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 140, 0, 6, 13, 14,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 149, 0, 6, 14, 15,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 158, 0, 6, 15, 16,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 167, 0, 6, 16, 17,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 176, 0, 6, 20, 21,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 185, 0, 6, 21, 22,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 194, 0, 6, 22, 23,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 203, 0, 6, 23, 24,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 212, 0, 6, 24, 25,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 221, 0, 6, 25, 26,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 230, 0, 6, 26, 27,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 239, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 242, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 245, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 248, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 251, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 254, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 257, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 260, 6, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 263, 6, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 266, 6, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 269, 6, 25, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 272, 6, 26, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 275, 6, 27, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 278, 6, 28, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 281, 6, 12, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 290, 6, 13, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 299, 6, 14, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 308, 6, 15, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 317, 6, 16, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 326, 6, 17, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 335, 6, 22, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 344, 6, 23, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 353, 6, 24, 65,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 362, 6, 25, 68,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 371, 6, 26, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 380, 6, 27, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 389, 6, 12, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 398, 6, 13, 80,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 407, 6, 14, 83,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 416, 6, 15, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 425, 6, 16, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 434, 6, 17, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 443, 6, 22, 95,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 452, 6, 23, 98,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 461, 6, 24, 101,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 470, 6, 25, 104,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 479, 6, 26, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 488, 6, 27, 110,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 497, 6, 35, 77,
                                                                       131, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 524, 6, 38, 80,
                                                                       140, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 551, 6, 41, 83,
                                                                       149, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 578, 6, 44, 86,
                                                                       158, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 605, 6, 47, 89,
                                                                       167, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 632, 6, 59, 95,
                                                                       194, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 659, 6, 62, 98,
                                                                       203, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 686, 6, 65, 101,
                                                                       212, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 713, 6, 68, 104,
                                                                       221, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 740, 6, 71, 107,
                                                                       230, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 767, 6, 10, 11,
                                                                       239, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 773, 6, 11, 12,
                                                                       242, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 779, 6, 12, 13,
                                                                       245, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 785, 6, 13, 14,
                                                                       248, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 791, 6, 14, 15,
                                                                       251, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 797, 6, 15, 16,
                                                                       254, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 803, 6, 16, 17,
                                                                       257, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 809, 6, 20, 21,
                                                                       260, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 815, 6, 21, 22,
                                                                       263, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 821, 6, 22, 23,
                                                                       266, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 827, 6, 23, 24,
                                                                       269, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 833, 6, 24, 25,
                                                                       272, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 839, 6, 25, 26,
                                                                       275, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 845, 6, 26, 27,
                                                                       278, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 851, 3, 6, 767,
                                                                       239, 773, 281, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 869, 3, 6, 773,
                                                                       242, 779, 290, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 887, 3, 6, 779,
                                                                       245, 785, 299, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 905, 3, 6, 785,
                                                                       248, 791, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 923, 3, 6, 791,
                                                                       251, 797, 317, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 941, 3, 6, 797,
                                                                       254, 803, 326, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 959, 3, 6, 809,
                                                                       260, 815, 335, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 977, 3, 6, 815,
                                                                       263, 821, 344, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 995, 3, 6, 821,
                                                                       266, 827, 353, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1013, 3, 6, 827,
                                                                       269, 833, 362, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1031, 3, 6, 833,
                                                                       272, 839, 371, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1049, 3, 6, 839,
                                                                       275, 845, 380, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1067, 0, 6, 767,
                                                                       239, 773, 389, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1085, 0, 6, 773,
                                                                       242, 779, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1103, 0, 6, 779,
                                                                       245, 785, 407, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1121, 0, 6, 785,
                                                                       248, 791, 416, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1139, 0, 6, 791,
                                                                       251, 797, 425, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1157, 0, 6, 797,
                                                                       254, 803, 434, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1175, 0, 6, 809,
                                                                       260, 815, 443, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1193, 0, 6, 815,
                                                                       263, 821, 452, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1211, 0, 6, 821,
                                                                       266, 827, 461, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1229, 0, 6, 827,
                                                                       269, 833, 470, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1247, 0, 6, 833,
                                                                       272, 839, 479, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1265, 0, 6, 839,
                                                                       275, 845, 488, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1283, 0, 3, 6,
                                                                       851, 281, 869, 1067, 389,
                                                                       1085, 113, 122, 497,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1337, 0, 3, 6,
                                                                       869, 290, 887, 1085, 398,
                                                                       1103, 122, 131, 524,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1391, 0, 3, 6,
                                                                       887, 299, 905, 1103, 407,
                                                                       1121, 131, 140, 551,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1445, 0, 3, 6,
                                                                       905, 308, 923, 1121, 416,
                                                                       1139, 140, 149, 578,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1499, 0, 3, 6,
                                                                       923, 317, 941, 1139, 425,
                                                                       1157, 149, 158, 605,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1553, 0, 3, 6,
                                                                       959, 335, 977, 1175, 443,
                                                                       1193, 176, 185, 632,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1607, 0, 3, 6,
                                                                       977, 344, 995, 1193, 452,
                                                                       1211, 185, 194, 659,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1661, 0, 3, 6,
                                                                       995, 353, 1013, 1211, 461,
                                                                       1229, 194, 203, 686,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1715, 0, 3, 6,
                                                                       1013, 362, 1031, 1229,
                                                                       470, 1247, 203, 212, 713,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1769, 0, 3, 6,
                                                                       1031, 371, 1049, 1247,
                                                                       479, 1265, 212, 221, 740,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1823, 6, 239, 242,
                                                                       779, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1833, 6, 242, 245,
                                                                       785, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1843, 6, 245, 248,
                                                                       791, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1853, 6, 248, 251,
                                                                       797, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1863, 6, 251, 254,
                                                                       803, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1873, 6, 260, 263,
                                                                       821, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1883, 6, 263, 266,
                                                                       827, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1893, 6, 266, 269,
                                                                       833, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1903, 6, 269, 272,
                                                                       839, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1913, 6, 272, 275,
                                                                       845, ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1923, 3, 6, 1823,
                                                                       779, 1833, 887, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1953, 3, 6, 1833,
                                                                       785, 1843, 905, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1983, 3, 6, 1843,
                                                                       791, 1853, 923, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2013, 3, 6, 1853,
                                                                       797, 1863, 941, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2043, 3, 6, 1873,
                                                                       821, 1883, 995, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2073, 3, 6, 1883,
                                                                       827, 1893, 1013, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2103, 3, 6, 1893,
                                                                       833, 1903, 1031, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 2133, 3, 6, 1903,
                                                                       839, 1913, 1049, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2163, 0, 6, 1823,
                                                                       779, 1833, 389, 398, 1103,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2193, 0, 6, 1833,
                                                                       785, 1843, 398, 407, 1121,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2223, 0, 6, 1843,
                                                                       791, 1853, 407, 416, 1139,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2253, 0, 6, 1853,
                                                                       797, 1863, 416, 425, 1157,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2283, 0, 6, 1873,
                                                                       821, 1883, 443, 452, 1211,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2313, 0, 6, 1883,
                                                                       827, 1893, 452, 461, 1229,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2343, 0, 6, 1893,
                                                                       833, 1903, 461, 470, 1247,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2373, 0, 6, 1903,
                                                                       839, 1913, 470, 479, 1265,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2403, 0, 3, 6,
                                                                       1923, 887, 1953, 2163,
                                                                       1103, 2193, 497, 524,
                                                                       1391, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2493, 0, 3, 6,
                                                                       1953, 905, 1983, 2193,
                                                                       1121, 2223, 524, 551,
                                                                       1445, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2583, 0, 3, 6,
                                                                       1983, 923, 2013, 2223,
                                                                       1139, 2253, 551, 578,
                                                                       1499, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2673, 0, 3, 6,
                                                                       2043, 995, 2073, 2283,
                                                                       1211, 2313, 632, 659,
                                                                       1661, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2763, 0, 3, 6,
                                                                       2073, 1013, 2103, 2313,
                                                                       1229, 2343, 659, 686,
                                                                       1715, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2853, 0, 3, 6,
                                                                       2103, 1031, 2133, 2343,
                                                                       1247, 2373, 686, 713,
                                                                       1769, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2943, 6, 767, 773,
                                                                       1823, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2958, 6, 773, 779,
                                                                       1833, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2973, 6, 779, 785,
                                                                       1843, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2988, 6, 785, 791,
                                                                       1853, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3003, 6, 791, 797,
                                                                       1863, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3018, 6, 809, 815,
                                                                       1873, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3033, 6, 815, 821,
                                                                       1883, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3048, 6, 821, 827,
                                                                       1893, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3063, 6, 827, 833,
                                                                       1903, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3078, 6, 833, 839,
                                                                       1913, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3093, 3, 6, 2943,
                                                                       1823, 2958, 851, 869,
                                                                       1923, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3138, 3, 6, 2958,
                                                                       1833, 2973, 869, 887,
                                                                       1953, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3183, 3, 6, 2973,
                                                                       1843, 2988, 887, 905,
                                                                       1983, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3228, 3, 6, 2988,
                                                                       1853, 3003, 905, 923,
                                                                       2013, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3273, 3, 6, 3018,
                                                                       1873, 3033, 959, 977,
                                                                       2043, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3318, 3, 6, 3033,
                                                                       1883, 3048, 977, 995,
                                                                       2073, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3363, 3, 6, 3048,
                                                                       1893, 3063, 995, 1013,
                                                                       2103, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 3408, 3, 6, 3063,
                                                                       1903, 3078, 1013, 1031,
                                                                       2133, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3453, 0, 6, 2943,
                                                                       1823, 2958, 1067, 1085,
                                                                       2163, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3498, 0, 6, 2958,
                                                                       1833, 2973, 1085, 1103,
                                                                       2193, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3543, 0, 6, 2973,
                                                                       1843, 2988, 1103, 1121,
                                                                       2223, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3588, 0, 6, 2988,
                                                                       1853, 3003, 1121, 1139,
                                                                       2253, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3633, 0, 6, 3018,
                                                                       1873, 3033, 1175, 1193,
                                                                       2283, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3678, 0, 6, 3033,
                                                                       1883, 3048, 1193, 1211,
                                                                       2313, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3723, 0, 6, 3048,
                                                                       1893, 3063, 1211, 1229,
                                                                       2343, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3768, 0, 6, 3063,
                                                                       1903, 3078, 1229, 1247,
                                                                       2373, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 3813, 0, 3, 6,
                                                                       3093, 1923, 3138, 3453,
                                                                       2163, 3498, 1283, 1337,
                                                                       2403, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 3948, 0, 3, 6,
                                                                       3138, 1953, 3183, 3498,
                                                                       2193, 3543, 1337, 1391,
                                                                       2493, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 4083, 0, 3, 6,
                                                                       3183, 1983, 3228, 3543,
                                                                       2223, 3588, 1391, 1445,
                                                                       2583, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 4218, 0, 3, 6,
                                                                       3273, 2043, 3318, 3633,
                                                                       2283, 3678, 1553, 1607,
                                                                       2673, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 4353, 0, 3, 6,
                                                                       3318, 2073, 3363, 3678,
                                                                       2313, 3723, 1607, 1661,
                                                                       2763, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 4488, 0, 3, 6,
                                                                       3363, 2103, 3408, 3723,
                                                                       2343, 3768, 1661, 1715,
                                                                       2853, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4623, 6, 1823,
                                                                       1833, 2973, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4644, 6, 1833,
                                                                       1843, 2988, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4665, 6, 1843,
                                                                       1853, 3003, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4686, 6, 1873,
                                                                       1883, 3048, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4707, 6, 1883,
                                                                       1893, 3063, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4728, 6, 1893,
                                                                       1903, 3078, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 4749, 3, 6, 4623,
                                                                       2973, 4644, 1923, 1953,
                                                                       3183, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 4812, 3, 6, 4644,
                                                                       2988, 4665, 1953, 1983,
                                                                       3228, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 4875, 3, 6, 4686,
                                                                       3048, 4707, 2043, 2073,
                                                                       3363, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 4938, 3, 6, 4707,
                                                                       3063, 4728, 2073, 2103,
                                                                       3408, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 5001, 0, 6, 4623,
                                                                       2973, 4644, 2163, 2193,
                                                                       3543, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 5064, 0, 6, 4644,
                                                                       2988, 4665, 2193, 2223,
                                                                       3588, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 5127, 0, 6, 4686,
                                                                       3048, 4707, 2283, 2313,
                                                                       3723, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 5190, 0, 6, 4707,
                                                                       3063, 4728, 2313, 2343,
                                                                       3768, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 5253, 0, 3, 6,
                                                                       4749, 3183, 4812, 5001,
                                                                       3543, 5064, 2403, 2493,
                                                                       4083, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 5442, 0, 3, 6,
                                                                       4875, 3363, 4938, 5127,
                                                                       3723, 5190, 2673, 2763,
                                                                       4488, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5631, 6, 2943,
                                                                       2958, 4623, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5659, 6, 2958,
                                                                       2973, 4644, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5687, 6, 2973,
                                                                       2988, 4665, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5715, 6, 3018,
                                                                       3033, 4686, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5743, 6, 3033,
                                                                       3048, 4707, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 5771, 6, 3048,
                                                                       3063, 4728, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 5799, 3, 6, 5631,
                                                                       4623, 5659, 3093, 3138,
                                                                       4749, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 5883, 3, 6, 5659,
                                                                       4644, 5687, 3138, 3183,
                                                                       4812, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 5967, 3, 6, 5715,
                                                                       4686, 5743, 3273, 3318,
                                                                       4875, ncols, gamma, p,
                                                                       q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 6051, 3, 6, 5743,
                                                                       4707, 5771, 3318, 3363,
                                                                       4938, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 6135, 0, 6, 5631,
                                                                       4623, 5659, 3453, 3498,
                                                                       5001, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 6219, 0, 6, 5659,
                                                                       4644, 5687, 3498, 3543,
                                                                       5064, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 6303, 0, 6, 5715,
                                                                       4686, 5743, 3633, 3678,
                                                                       5127, ncols, gamma, p,
                                                                       q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 6387, 0, 6, 5743,
                                                                       4707, 5771, 3678, 3723,
                                                                       5190, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 6471, 0, 3, 6,
                                                                       5799, 4749, 5883, 6135,
                                                                       5001, 6219, 3813, 3948,
                                                                       5253, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 6723, 0, 3, 6,
                                                                       5967, 4875, 6051, 6303,
                                                                       5127, 6387, 4218, 4353,
                                                                       5442, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_s_x(buffer, 6975, 6723, 3, 28, ncols, beta);

                    simdgeo::geom_s_y(buffer, 7059, 6723, 3, 28, ncols, beta);

                    simdgeo::geom_s_z(buffer, 7143, 6723, 3, 28, ncols, beta);

                    simdgeo::geom_s_x(buffer, 7227, 6471, 3, 28, ncols, beta);

                    simdgeo::geom_s_y(buffer, 7311, 6471, 3, 28, ncols, beta);

                    simdgeo::geom_s_z(buffer, 7395, 6471, 3, 28, ncols, beta);

                    simdfunc::contract_primitives(buffer, 7479, 6975, 504, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 7983, 7479, 3, 1, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 7983, 13, nmax);

        simdtrf::transform_i_inner(buffer, 7983, 7563, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 39 * nvalues + n * npairs, nvalues, buffer, 7983, 13,
                                   nmax);

        simdtrf::transform_i_inner(buffer, 7983, 7647, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 78 * nvalues + n * npairs, nvalues, buffer, 7983, 13,
                                   nmax);

        simdtrf::transform_i_inner(buffer, 7983, 7731, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 117 * nvalues + n * npairs, nvalues, buffer, 7983,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 7983, 7815, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 156 * nvalues + n * npairs, nvalues, buffer, 7983,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 7983, 7899, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 195 * nvalues + n * npairs, nvalues, buffer, 7983,
                                   13, nmax);
    }

    for (size_t m = 0; m < 234; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
