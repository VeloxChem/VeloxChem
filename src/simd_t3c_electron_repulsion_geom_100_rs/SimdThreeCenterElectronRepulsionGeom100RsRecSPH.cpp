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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecSPH.hpp"

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
#include "SimdTransformH.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_sph_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_sph_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 4988, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 198 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 4988, 4577, 378, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto alpha = a_exps[i];

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

                    simdfunc::compute_t3c_erf_boys_function(buffer, coordinates, 9, 6, {1, 2, 3,
                                                            4, 5, 6, 7}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 17, 6, {1, 2, 3, 4,
                                                        5, 6, 7}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 3, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 3, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 3, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 3, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 3, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 67, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 70, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 73, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 76, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 79, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 82, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 85, 0, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 88, 0, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 91, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 94, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 97, 0, 6, 10, 11,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 106, 0, 6, 11, 12,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 115, 0, 6, 12, 13,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 124, 0, 6, 13, 14,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 133, 0, 6, 14, 15,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 142, 0, 6, 18, 19,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 151, 0, 6, 19, 20,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 160, 0, 6, 20, 21,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 169, 0, 6, 21, 22,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 178, 0, 6, 22, 23,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 187, 6, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 190, 6, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 193, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 196, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 199, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 202, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 205, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 208, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 211, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 214, 6, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 217, 6, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 220, 6, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 223, 6, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 226, 6, 24, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 229, 6, 12, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 238, 6, 13, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 247, 6, 14, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 256, 6, 15, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 265, 6, 20, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 274, 6, 21, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 283, 6, 22, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 292, 6, 23, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 301, 6, 10, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 310, 6, 11, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 319, 6, 12, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 328, 6, 13, 70,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 337, 6, 14, 73,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 346, 6, 15, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 355, 6, 18, 79,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 364, 6, 19, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 373, 6, 20, 85,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 382, 6, 21, 88,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 391, 6, 22, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 400, 6, 23, 94,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 409, 6, 25, 61,
                                                                       97, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 436, 6, 28, 64,
                                                                       106, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 463, 6, 31, 67,
                                                                       115, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 490, 6, 34, 70,
                                                                       124, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 517, 6, 37, 73,
                                                                       133, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 544, 6, 43, 79,
                                                                       142, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 571, 6, 46, 82,
                                                                       151, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 598, 6, 49, 85,
                                                                       160, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 625, 6, 52, 88,
                                                                       169, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 652, 6, 55, 91,
                                                                       178, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 679, 6, 10, 11,
                                                                       193, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 685, 6, 11, 12,
                                                                       196, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 691, 6, 12, 13,
                                                                       199, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 697, 6, 13, 14,
                                                                       202, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 703, 6, 14, 15,
                                                                       205, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 709, 6, 18, 19,
                                                                       214, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 715, 6, 19, 20,
                                                                       217, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 721, 6, 20, 21,
                                                                       220, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 727, 6, 21, 22,
                                                                       223, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 733, 6, 22, 23,
                                                                       226, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 739, 3, 6, 679,
                                                                       193, 685, 229, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 757, 3, 6, 685,
                                                                       196, 691, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 775, 3, 6, 691,
                                                                       199, 697, 247, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 793, 3, 6, 697,
                                                                       202, 703, 256, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 811, 3, 6, 709,
                                                                       214, 715, 265, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 829, 3, 6, 715,
                                                                       217, 721, 274, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 847, 3, 6, 721,
                                                                       220, 727, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 865, 3, 6, 727,
                                                                       223, 733, 292, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 883, 0, 6, 679,
                                                                       193, 685, 319, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 901, 0, 6, 685,
                                                                       196, 691, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 919, 0, 6, 691,
                                                                       199, 697, 337, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 937, 0, 6, 697,
                                                                       202, 703, 346, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 955, 0, 6, 709,
                                                                       214, 715, 373, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 973, 0, 6, 715,
                                                                       217, 721, 382, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 991, 0, 6, 721,
                                                                       220, 727, 391, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1009, 0, 6, 727,
                                                                       223, 733, 400, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1027, 0, 3, 6,
                                                                       739, 229, 757, 883, 319,
                                                                       901, 97, 106, 463, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1081, 0, 3, 6,
                                                                       757, 238, 775, 901, 328,
                                                                       919, 106, 115, 490, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1135, 0, 3, 6,
                                                                       775, 247, 793, 919, 337,
                                                                       937, 115, 124, 517, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1189, 0, 3, 6,
                                                                       811, 265, 829, 955, 373,
                                                                       973, 142, 151, 598, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1243, 0, 3, 6,
                                                                       829, 274, 847, 973, 382,
                                                                       991, 151, 160, 625, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1297, 0, 3, 6,
                                                                       847, 283, 865, 991, 391,
                                                                       1009, 160, 169, 652,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1351, 6, 187, 190,
                                                                       679, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1361, 6, 190, 193,
                                                                       685, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1371, 6, 193, 196,
                                                                       691, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1381, 6, 196, 199,
                                                                       697, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1391, 6, 199, 202,
                                                                       703, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1401, 6, 208, 211,
                                                                       709, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1411, 6, 211, 214,
                                                                       715, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1421, 6, 214, 217,
                                                                       721, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1431, 6, 217, 220,
                                                                       727, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1441, 6, 220, 223,
                                                                       733, ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1451, 3, 6, 1351,
                                                                       679, 1361, 739, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1481, 3, 6, 1361,
                                                                       685, 1371, 757, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1511, 3, 6, 1371,
                                                                       691, 1381, 775, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1541, 3, 6, 1381,
                                                                       697, 1391, 793, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1571, 3, 6, 1401,
                                                                       709, 1411, 811, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1601, 3, 6, 1411,
                                                                       715, 1421, 829, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1631, 3, 6, 1421,
                                                                       721, 1431, 847, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 1661, 3, 6, 1431,
                                                                       727, 1441, 865, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1691, 0, 6, 1351,
                                                                       679, 1361, 301, 310, 883,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1721, 0, 6, 1361,
                                                                       685, 1371, 310, 319, 901,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1751, 0, 6, 1371,
                                                                       691, 1381, 319, 328, 919,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1781, 0, 6, 1381,
                                                                       697, 1391, 328, 337, 937,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1811, 0, 6, 1401,
                                                                       709, 1411, 355, 364, 955,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1841, 0, 6, 1411,
                                                                       715, 1421, 364, 373, 973,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1871, 0, 6, 1421,
                                                                       721, 1431, 373, 382, 991,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1901, 0, 6, 1431,
                                                                       727, 1441, 382, 391, 1009,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 1931, 0, 3, 6,
                                                                       1451, 739, 1481, 1691,
                                                                       883, 1721, 409, 436, 1027,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2021, 0, 3, 6,
                                                                       1481, 757, 1511, 1721,
                                                                       901, 1751, 436, 463, 1081,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2111, 0, 3, 6,
                                                                       1511, 775, 1541, 1751,
                                                                       919, 1781, 463, 490, 1135,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2201, 0, 3, 6,
                                                                       1571, 811, 1601, 1811,
                                                                       955, 1841, 544, 571, 1189,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2291, 0, 3, 6,
                                                                       1601, 829, 1631, 1841,
                                                                       973, 1871, 571, 598, 1243,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 2381, 0, 3, 6,
                                                                       1631, 847, 1661, 1871,
                                                                       991, 1901, 598, 625, 1297,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2471, 6, 679, 685,
                                                                       1371, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2486, 6, 685, 691,
                                                                       1381, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2501, 6, 691, 697,
                                                                       1391, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2516, 6, 709, 715,
                                                                       1421, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2531, 6, 715, 721,
                                                                       1431, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 2546, 6, 721, 727,
                                                                       1441, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2561, 3, 6, 2471,
                                                                       1371, 2486, 739, 757,
                                                                       1511, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2606, 3, 6, 2486,
                                                                       1381, 2501, 757, 775,
                                                                       1541, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2651, 3, 6, 2516,
                                                                       1421, 2531, 811, 829,
                                                                       1631, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 2696, 3, 6, 2531,
                                                                       1431, 2546, 829, 847,
                                                                       1661, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 2741, 0, 6, 2471,
                                                                       1371, 2486, 883, 901,
                                                                       1751, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 2786, 0, 6, 2486,
                                                                       1381, 2501, 901, 919,
                                                                       1781, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 2831, 0, 6, 2516,
                                                                       1421, 2531, 955, 973,
                                                                       1871, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 2876, 0, 6, 2531,
                                                                       1431, 2546, 973, 991,
                                                                       1901, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 2921, 0, 3, 6,
                                                                       2561, 1511, 2606, 2741,
                                                                       1751, 2786, 1027, 1081,
                                                                       2111, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 3056, 0, 3, 6,
                                                                       2651, 1631, 2696, 2831,
                                                                       1871, 2876, 1189, 1243,
                                                                       2381, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3191, 6, 1351,
                                                                       1361, 2471, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3212, 6, 1361,
                                                                       1371, 2486, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3233, 6, 1371,
                                                                       1381, 2501, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3254, 6, 1401,
                                                                       1411, 2516, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3275, 6, 1411,
                                                                       1421, 2531, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 3296, 6, 1421,
                                                                       1431, 2546, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 3317, 3, 6, 3191,
                                                                       2471, 3212, 1451, 1481,
                                                                       2561, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 3380, 3, 6, 3212,
                                                                       2486, 3233, 1481, 1511,
                                                                       2606, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 3443, 3, 6, 3254,
                                                                       2516, 3275, 1571, 1601,
                                                                       2651, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 3506, 3, 6, 3275,
                                                                       2531, 3296, 1601, 1631,
                                                                       2696, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 3569, 0, 6, 3191,
                                                                       2471, 3212, 1691, 1721,
                                                                       2741, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 3632, 0, 6, 3212,
                                                                       2486, 3233, 1721, 1751,
                                                                       2786, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 3695, 0, 6, 3254,
                                                                       2516, 3275, 1811, 1841,
                                                                       2831, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 3758, 0, 6, 3275,
                                                                       2531, 3296, 1841, 1871,
                                                                       2876, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 3821, 0, 3, 6,
                                                                       3317, 2561, 3380, 3569,
                                                                       2741, 3632, 1931, 2021,
                                                                       2921, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 4010, 0, 3, 6,
                                                                       3443, 2651, 3506, 3695,
                                                                       2831, 3758, 2201, 2291,
                                                                       3056, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_s_x(buffer, 4199, 4010, 1, 63, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 4262, 4010, 1, 63, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 4325, 4010, 1, 63, ncols, alpha);

                    simdgeo::geom_s_x(buffer, 4388, 3821, 1, 63, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 4451, 3821, 1, 63, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 4514, 3821, 1, 63, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 4577, 4199, 378, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 4955, 4577, 3, 1, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 4955, 11, nmax);

        simdtrf::transform_h_inner(buffer, 4955, 4640, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 33 * nvalues + n * npairs, nvalues, buffer, 4955, 11,
                                   nmax);

        simdtrf::transform_h_inner(buffer, 4955, 4703, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 66 * nvalues + n * npairs, nvalues, buffer, 4955, 11,
                                   nmax);

        simdtrf::transform_h_inner(buffer, 4955, 4766, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 99 * nvalues + n * npairs, nvalues, buffer, 4955, 11,
                                   nmax);

        simdtrf::transform_h_inner(buffer, 4955, 4829, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 132 * nvalues + n * npairs, nvalues, buffer, 4955,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 4955, 4892, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 165 * nvalues + n * npairs, nvalues, buffer, 4955,
                                   11, nmax);
    }

    for (size_t m = 0; m < 198; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
