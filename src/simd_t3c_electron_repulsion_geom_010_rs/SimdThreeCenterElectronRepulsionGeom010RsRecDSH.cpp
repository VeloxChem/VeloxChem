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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecDSH.hpp"

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
#include "SimdTransformD.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_dsh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_dsh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 13233, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 330 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 13233, 12411, 756, dimensions);

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

                    simdfunc::compute_t3c_erf_boys_function(buffer, coordinates, 9, 6, {1, 2, 3,
                                                            4, 5, 6, 7, 8}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 18, 6, {1, 2, 3, 4,
                                                        5, 6, 7, 8}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 3, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 3, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 3, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 3, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 60, 3, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 63, 3, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 66, 3, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 69, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 72, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 75, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 78, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 81, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 84, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 87, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 90, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 93, 0, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 96, 0, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 99, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 102, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 105, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 108, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 111, 0, 6, 10, 11,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 120, 0, 6, 11, 12,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 129, 0, 6, 12, 13,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 138, 0, 6, 13, 14,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 147, 0, 6, 14, 15,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 156, 0, 6, 15, 16,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 165, 0, 6, 19, 20,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 174, 0, 6, 20, 21,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 183, 0, 6, 21, 22,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 192, 0, 6, 22, 23,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 201, 0, 6, 23, 24,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 210, 0, 6, 24, 25,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 219, 0, 6, 10, 11,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 225, 0, 6, 11, 12,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 231, 0, 6, 12, 13,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 237, 0, 6, 13, 14,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 243, 0, 6, 14, 15,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 249, 0, 6, 15, 16,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 255, 0, 6, 19, 20,
                                                                       90, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 261, 0, 6, 20, 21,
                                                                       93, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 267, 0, 6, 21, 22,
                                                                       96, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 273, 0, 6, 22, 23,
                                                                       99, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 279, 0, 6, 23, 24,
                                                                       102, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 285, 0, 6, 24, 25,
                                                                       105, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 291, 0, 3, 6, 69,
                                                                       72, 111, 120, 219, 225,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 309, 0, 3, 6, 72,
                                                                       75, 120, 129, 225, 231,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 327, 0, 3, 6, 75,
                                                                       78, 129, 138, 231, 237,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 345, 0, 3, 6, 78,
                                                                       81, 138, 147, 237, 243,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 363, 0, 3, 6, 81,
                                                                       84, 147, 156, 243, 249,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 381, 0, 3, 6, 90,
                                                                       93, 165, 174, 255, 261,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 399, 0, 3, 6, 93,
                                                                       96, 174, 183, 261, 267,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 417, 0, 3, 6, 96,
                                                                       99, 183, 192, 267, 273,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 435, 0, 3, 6, 99,
                                                                       102, 192, 201, 273, 279,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 453, 0, 3, 6, 102,
                                                                       105, 201, 210, 279, 285,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 471, 6, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 474, 6, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 477, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 480, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 483, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 486, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 489, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 492, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 495, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 498, 6, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 501, 6, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 504, 6, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 507, 6, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 510, 6, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 513, 6, 25, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 516, 6, 26, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 519, 6, 12, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 528, 6, 13, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 537, 6, 14, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 546, 6, 15, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 555, 6, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 564, 6, 21, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 573, 6, 22, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 582, 6, 23, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 591, 6, 24, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 600, 6, 25, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 609, 6, 10, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 618, 6, 11, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 627, 6, 12, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 636, 6, 13, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 645, 6, 14, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 654, 6, 15, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 663, 6, 16, 87,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 672, 6, 19, 90,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 681, 6, 20, 93,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 690, 6, 21, 96,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 699, 6, 22, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 708, 6, 23, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 717, 6, 24, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 726, 6, 25, 108,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 735, 6, 27, 69,
                                                                       111, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 762, 6, 30, 72,
                                                                       120, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 789, 6, 33, 75,
                                                                       129, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 816, 6, 36, 78,
                                                                       138, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 843, 6, 39, 81,
                                                                       147, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 870, 6, 42, 84,
                                                                       156, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 897, 6, 48, 90,
                                                                       165, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 924, 6, 51, 93,
                                                                       174, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 951, 6, 54, 96,
                                                                       183, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 978, 6, 57, 99,
                                                                       192, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1005, 6, 60, 102,
                                                                       201, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1032, 6, 63, 105,
                                                                       210, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1059, 6, 69, 219,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1077, 6, 72, 225,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1095, 6, 75, 231,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1113, 6, 78, 237,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1131, 6, 81, 243,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1149, 6, 84, 249,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1167, 6, 90, 255,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1185, 6, 93, 261,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1203, 6, 96, 267,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1221, 6, 99, 273,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1239, 6, 102, 279,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1257, 6, 105, 285,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1275, 0, 6, 735,
                                                                       111, 762, 219, 291, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1329, 0, 6, 762,
                                                                       120, 789, 225, 309, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1383, 0, 6, 789,
                                                                       129, 816, 231, 327, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1437, 0, 6, 816,
                                                                       138, 843, 237, 345, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1491, 0, 6, 843,
                                                                       147, 870, 243, 363, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1545, 0, 6, 897,
                                                                       165, 924, 255, 381, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1599, 0, 6, 924,
                                                                       174, 951, 261, 399, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1653, 0, 6, 951,
                                                                       183, 978, 267, 417, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1707, 0, 6, 978,
                                                                       192, 1005, 273, 435,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1761, 0, 6, 1005,
                                                                       201, 1032, 279, 453,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1815, 6, 10, 11,
                                                                       477, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1821, 6, 11, 12,
                                                                       480, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1827, 6, 12, 13,
                                                                       483, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1833, 6, 13, 14,
                                                                       486, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1839, 6, 14, 15,
                                                                       489, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1845, 6, 15, 16,
                                                                       492, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1851, 6, 19, 20,
                                                                       501, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1857, 6, 20, 21,
                                                                       504, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1863, 6, 21, 22,
                                                                       507, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1869, 6, 22, 23,
                                                                       510, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1875, 6, 23, 24,
                                                                       513, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1881, 6, 24, 25,
                                                                       516, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1887, 3, 6, 1815,
                                                                       477, 1821, 519, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1905, 3, 6, 1821,
                                                                       480, 1827, 528, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1923, 3, 6, 1827,
                                                                       483, 1833, 537, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1941, 3, 6, 1833,
                                                                       486, 1839, 546, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1959, 3, 6, 1839,
                                                                       489, 1845, 555, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1977, 3, 6, 1851,
                                                                       501, 1857, 564, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1995, 3, 6, 1857,
                                                                       504, 1863, 573, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2013, 3, 6, 1863,
                                                                       507, 1869, 582, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2031, 3, 6, 1869,
                                                                       510, 1875, 591, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2049, 3, 6, 1875,
                                                                       513, 1881, 600, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2067, 0, 6, 1815,
                                                                       477, 1821, 627, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2085, 0, 6, 1821,
                                                                       480, 1827, 636, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2103, 0, 6, 1827,
                                                                       483, 1833, 645, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2121, 0, 6, 1833,
                                                                       486, 1839, 654, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2139, 0, 6, 1839,
                                                                       489, 1845, 663, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2157, 0, 6, 1851,
                                                                       501, 1857, 690, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2175, 0, 6, 1857,
                                                                       504, 1863, 699, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2193, 0, 6, 1863,
                                                                       507, 1869, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2211, 0, 6, 1869,
                                                                       510, 1875, 717, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2229, 0, 6, 1875,
                                                                       513, 1881, 726, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2247, 0, 3, 6,
                                                                       1887, 519, 1905, 2067,
                                                                       627, 2085, 111, 120, 789,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2301, 0, 3, 6,
                                                                       1905, 528, 1923, 2085,
                                                                       636, 2103, 120, 129, 816,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2355, 0, 3, 6,
                                                                       1923, 537, 1941, 2103,
                                                                       645, 2121, 129, 138, 843,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2409, 0, 3, 6,
                                                                       1941, 546, 1959, 2121,
                                                                       654, 2139, 138, 147, 870,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2463, 0, 3, 6,
                                                                       1977, 564, 1995, 2157,
                                                                       690, 2175, 165, 174, 951,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2517, 0, 3, 6,
                                                                       1995, 573, 2013, 2175,
                                                                       699, 2193, 174, 183, 978,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2571, 0, 3, 6,
                                                                       2013, 582, 2031, 2193,
                                                                       708, 2211, 183, 192, 1005,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2625, 0, 3, 6,
                                                                       2031, 591, 2049, 2211,
                                                                       717, 2229, 192, 201, 1032,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2679, 0, 6, 2067,
                                                                       627, 2085, 219, 225, 1095,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2715, 0, 6, 2085,
                                                                       636, 2103, 225, 231, 1113,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2751, 0, 6, 2103,
                                                                       645, 2121, 231, 237, 1131,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2787, 0, 6, 2121,
                                                                       654, 2139, 237, 243, 1149,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2823, 0, 6, 2157,
                                                                       690, 2175, 255, 261, 1203,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2859, 0, 6, 2175,
                                                                       699, 2193, 261, 267, 1221,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2895, 0, 6, 2193,
                                                                       708, 2211, 267, 273, 1239,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2931, 0, 6, 2211,
                                                                       717, 2229, 273, 279, 1257,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 2967, 0, 3, 6,
                                                                       2247, 789, 2301, 2679,
                                                                       1095, 2715, 291, 309,
                                                                       1383, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3075, 0, 3, 6,
                                                                       2301, 816, 2355, 2715,
                                                                       1113, 2751, 309, 327,
                                                                       1437, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3183, 0, 3, 6,
                                                                       2355, 843, 2409, 2751,
                                                                       1131, 2787, 327, 345,
                                                                       1491, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3291, 0, 3, 6,
                                                                       2463, 951, 2517, 2823,
                                                                       1203, 2859, 381, 399,
                                                                       1653, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3399, 0, 3, 6,
                                                                       2517, 978, 2571, 2859,
                                                                       1221, 2895, 399, 417,
                                                                       1707, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3507, 0, 3, 6,
                                                                       2571, 1005, 2625, 2895,
                                                                       1239, 2931, 417, 435,
                                                                       1761, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3615, 6, 471, 474,
                                                                       1815, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3625, 6, 474, 477,
                                                                       1821, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3635, 6, 477, 480,
                                                                       1827, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3645, 6, 480, 483,
                                                                       1833, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3655, 6, 483, 486,
                                                                       1839, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3665, 6, 486, 489,
                                                                       1845, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3675, 6, 495, 498,
                                                                       1851, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3685, 6, 498, 501,
                                                                       1857, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3695, 6, 501, 504,
                                                                       1863, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3705, 6, 504, 507,
                                                                       1869, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3715, 6, 507, 510,
                                                                       1875, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3725, 6, 510, 513,
                                                                       1881, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3735, 3, 6, 3615,
                                                                       1815, 3625, 1887, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3765, 3, 6, 3625,
                                                                       1821, 3635, 1905, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3795, 3, 6, 3635,
                                                                       1827, 3645, 1923, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3825, 3, 6, 3645,
                                                                       1833, 3655, 1941, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3855, 3, 6, 3655,
                                                                       1839, 3665, 1959, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3885, 3, 6, 3675,
                                                                       1851, 3685, 1977, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3915, 3, 6, 3685,
                                                                       1857, 3695, 1995, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3945, 3, 6, 3695,
                                                                       1863, 3705, 2013, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3975, 3, 6, 3705,
                                                                       1869, 3715, 2031, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4005, 3, 6, 3715,
                                                                       1875, 3725, 2049, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4035, 0, 6, 3615,
                                                                       1815, 3625, 609, 618,
                                                                       2067, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4065, 0, 6, 3625,
                                                                       1821, 3635, 618, 627,
                                                                       2085, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4095, 0, 6, 3635,
                                                                       1827, 3645, 627, 636,
                                                                       2103, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4125, 0, 6, 3645,
                                                                       1833, 3655, 636, 645,
                                                                       2121, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4155, 0, 6, 3655,
                                                                       1839, 3665, 645, 654,
                                                                       2139, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4185, 0, 6, 3675,
                                                                       1851, 3685, 672, 681,
                                                                       2157, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4215, 0, 6, 3685,
                                                                       1857, 3695, 681, 690,
                                                                       2175, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4245, 0, 6, 3695,
                                                                       1863, 3705, 690, 699,
                                                                       2193, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4275, 0, 6, 3705,
                                                                       1869, 3715, 699, 708,
                                                                       2211, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4305, 0, 6, 3715,
                                                                       1875, 3725, 708, 717,
                                                                       2229, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4335, 0, 3, 6,
                                                                       3735, 1887, 3765, 4035,
                                                                       2067, 4065, 735, 762,
                                                                       2247, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4425, 0, 3, 6,
                                                                       3765, 1905, 3795, 4065,
                                                                       2085, 4095, 762, 789,
                                                                       2301, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4515, 0, 3, 6,
                                                                       3795, 1923, 3825, 4095,
                                                                       2103, 4125, 789, 816,
                                                                       2355, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4605, 0, 3, 6,
                                                                       3825, 1941, 3855, 4125,
                                                                       2121, 4155, 816, 843,
                                                                       2409, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4695, 0, 3, 6,
                                                                       3885, 1977, 3915, 4185,
                                                                       2157, 4215, 897, 924,
                                                                       2463, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4785, 0, 3, 6,
                                                                       3915, 1995, 3945, 4215,
                                                                       2175, 4245, 924, 951,
                                                                       2517, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4875, 0, 3, 6,
                                                                       3945, 2013, 3975, 4245,
                                                                       2193, 4275, 951, 978,
                                                                       2571, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4965, 0, 3, 6,
                                                                       3975, 2031, 4005, 4275,
                                                                       2211, 4305, 978, 1005,
                                                                       2625, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5055, 0, 6, 4035,
                                                                       2067, 4065, 1059, 1077,
                                                                       2679, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5115, 0, 6, 4065,
                                                                       2085, 4095, 1077, 1095,
                                                                       2715, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5175, 0, 6, 4095,
                                                                       2103, 4125, 1095, 1113,
                                                                       2751, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5235, 0, 6, 4125,
                                                                       2121, 4155, 1113, 1131,
                                                                       2787, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5295, 0, 6, 4185,
                                                                       2157, 4215, 1167, 1185,
                                                                       2823, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5355, 0, 6, 4215,
                                                                       2175, 4245, 1185, 1203,
                                                                       2859, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5415, 0, 6, 4245,
                                                                       2193, 4275, 1203, 1221,
                                                                       2895, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5475, 0, 6, 4275,
                                                                       2211, 4305, 1221, 1239,
                                                                       2931, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 5535, 0, 3, 6,
                                                                       4335, 2247, 4425, 5055,
                                                                       2679, 5115, 1275, 1329,
                                                                       2967, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 5715, 0, 3, 6,
                                                                       4425, 2301, 4515, 5115,
                                                                       2715, 5175, 1329, 1383,
                                                                       3075, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 5895, 0, 3, 6,
                                                                       4515, 2355, 4605, 5175,
                                                                       2751, 5235, 1383, 1437,
                                                                       3183, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 6075, 0, 3, 6,
                                                                       4695, 2463, 4785, 5295,
                                                                       2823, 5355, 1545, 1599,
                                                                       3291, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 6255, 0, 3, 6,
                                                                       4785, 2517, 4875, 5355,
                                                                       2859, 5415, 1599, 1653,
                                                                       3399, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 6435, 0, 3, 6,
                                                                       4875, 2571, 4965, 5415,
                                                                       2895, 5475, 1653, 1707,
                                                                       3507, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6615, 6, 1815,
                                                                       1821, 3635, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6630, 6, 1821,
                                                                       1827, 3645, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6645, 6, 1827,
                                                                       1833, 3655, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6660, 6, 1833,
                                                                       1839, 3665, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6675, 6, 1851,
                                                                       1857, 3695, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6690, 6, 1857,
                                                                       1863, 3705, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6705, 6, 1863,
                                                                       1869, 3715, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 6720, 6, 1869,
                                                                       1875, 3725, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6735, 3, 6, 6615,
                                                                       3635, 6630, 1887, 1905,
                                                                       3795, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6780, 3, 6, 6630,
                                                                       3645, 6645, 1905, 1923,
                                                                       3825, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6825, 3, 6, 6645,
                                                                       3655, 6660, 1923, 1941,
                                                                       3855, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6870, 3, 6, 6675,
                                                                       3695, 6690, 1977, 1995,
                                                                       3945, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6915, 3, 6, 6690,
                                                                       3705, 6705, 1995, 2013,
                                                                       3975, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 6960, 3, 6, 6705,
                                                                       3715, 6720, 2013, 2031,
                                                                       4005, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7005, 0, 6, 6615,
                                                                       3635, 6630, 2067, 2085,
                                                                       4095, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7050, 0, 6, 6630,
                                                                       3645, 6645, 2085, 2103,
                                                                       4125, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7095, 0, 6, 6645,
                                                                       3655, 6660, 2103, 2121,
                                                                       4155, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7140, 0, 6, 6675,
                                                                       3695, 6690, 2157, 2175,
                                                                       4245, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7185, 0, 6, 6690,
                                                                       3705, 6705, 2175, 2193,
                                                                       4275, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7230, 0, 6, 6705,
                                                                       3715, 6720, 2193, 2211,
                                                                       4305, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 7275, 0, 3, 6,
                                                                       6735, 3795, 6780, 7005,
                                                                       4095, 7050, 2247, 2301,
                                                                       4515, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 7410, 0, 3, 6,
                                                                       6780, 3825, 6825, 7050,
                                                                       4125, 7095, 2301, 2355,
                                                                       4605, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 7545, 0, 3, 6,
                                                                       6870, 3945, 6915, 7140,
                                                                       4245, 7185, 2463, 2517,
                                                                       4875, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 7680, 0, 3, 6,
                                                                       6915, 3975, 6960, 7185,
                                                                       4275, 7230, 2517, 2571,
                                                                       4965, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 7815, 0, 6, 7005,
                                                                       4095, 7050, 2679, 2715,
                                                                       5175, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 7905, 0, 6, 7050,
                                                                       4125, 7095, 2715, 2751,
                                                                       5235, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 7995, 0, 6, 7140,
                                                                       4245, 7185, 2823, 2859,
                                                                       5415, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8085, 0, 6, 7185,
                                                                       4275, 7230, 2859, 2895,
                                                                       5475, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 8175, 0, 3, 6,
                                                                       7275, 4515, 7410, 7815,
                                                                       5175, 7905, 2967, 3075,
                                                                       5895, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 8445, 0, 3, 6,
                                                                       7545, 4875, 7680, 7995,
                                                                       5415, 8085, 3291, 3399,
                                                                       6435, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8715, 6, 3615,
                                                                       3625, 6615, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8736, 6, 3625,
                                                                       3635, 6630, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8757, 6, 3635,
                                                                       3645, 6645, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8778, 6, 3645,
                                                                       3655, 6660, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8799, 6, 3675,
                                                                       3685, 6675, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8820, 6, 3685,
                                                                       3695, 6690, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8841, 6, 3695,
                                                                       3705, 6705, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 8862, 6, 3705,
                                                                       3715, 6720, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 8883, 3, 6, 8715,
                                                                       6615, 8736, 3735, 3765,
                                                                       6735, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 8946, 3, 6, 8736,
                                                                       6630, 8757, 3765, 3795,
                                                                       6780, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9009, 3, 6, 8757,
                                                                       6645, 8778, 3795, 3825,
                                                                       6825, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9072, 3, 6, 8799,
                                                                       6675, 8820, 3885, 3915,
                                                                       6870, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9135, 3, 6, 8820,
                                                                       6690, 8841, 3915, 3945,
                                                                       6915, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9198, 3, 6, 8841,
                                                                       6705, 8862, 3945, 3975,
                                                                       6960, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9261, 0, 6, 8715,
                                                                       6615, 8736, 4035, 4065,
                                                                       7005, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9324, 0, 6, 8736,
                                                                       6630, 8757, 4065, 4095,
                                                                       7050, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9387, 0, 6, 8757,
                                                                       6645, 8778, 4095, 4125,
                                                                       7095, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9450, 0, 6, 8799,
                                                                       6675, 8820, 4185, 4215,
                                                                       7140, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9513, 0, 6, 8820,
                                                                       6690, 8841, 4215, 4245,
                                                                       7185, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9576, 0, 6, 8841,
                                                                       6705, 8862, 4245, 4275,
                                                                       7230, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 9639, 0, 3, 6,
                                                                       8883, 6735, 8946, 9261,
                                                                       7005, 9324, 4335, 4425,
                                                                       7275, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 9828, 0, 3, 6,
                                                                       8946, 6780, 9009, 9324,
                                                                       7050, 9387, 4425, 4515,
                                                                       7410, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 10017, 0, 3, 6,
                                                                       9072, 6870, 9135, 9450,
                                                                       7140, 9513, 4695, 4785,
                                                                       7545, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 10206, 0, 3, 6,
                                                                       9135, 6915, 9198, 9513,
                                                                       7185, 9576, 4785, 4875,
                                                                       7680, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 10395, 0, 6, 9261,
                                                                       7005, 9324, 5055, 5115,
                                                                       7815, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 10521, 0, 6, 9324,
                                                                       7050, 9387, 5115, 5175,
                                                                       7905, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 10647, 0, 6, 9450,
                                                                       7140, 9513, 5295, 5355,
                                                                       7995, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 10773, 0, 6, 9513,
                                                                       7185, 9576, 5355, 5415,
                                                                       8085, ncols, gamma, p,
                                                                       q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 10899, 0, 3, 6,
                                                                       9639, 7275, 9828, 10395,
                                                                       7815, 10521, 5535, 5715,
                                                                       8175, ncols, gamma, p,
                                                                       q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 11277, 0, 3, 6,
                                                                       10017, 7545, 10206, 10647,
                                                                       7995, 10773, 6075, 6255,
                                                                       8445, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_s_x(buffer, 11655, 11277, 6, 21, ncols, beta);

                    simdgeo::geom_s_y(buffer, 11781, 11277, 6, 21, ncols, beta);

                    simdgeo::geom_s_z(buffer, 11907, 11277, 6, 21, ncols, beta);

                    simdgeo::geom_s_x(buffer, 12033, 10899, 6, 21, ncols, beta);

                    simdgeo::geom_s_y(buffer, 12159, 10899, 6, 21, ncols, beta);

                    simdgeo::geom_s_z(buffer, 12285, 10899, 6, 21, ncols, beta);

                    simdfunc::contract_primitives(buffer, 12411, 11655, 756, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 13167, 12411, 6, 1, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 13167, 11, nmax);

        simdtrf::transform_h_inner(buffer, 13167, 12537, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 55 * nvalues + n * npairs, nvalues, buffer, 13167,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 13167, 12663, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 110 * nvalues + n * npairs, nvalues, buffer, 13167,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 13167, 12789, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 165 * nvalues + n * npairs, nvalues, buffer, 13167,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 13167, 12915, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 220 * nvalues + n * npairs, nvalues, buffer, 13167,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 13167, 13041, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 275 * nvalues + n * npairs, nvalues, buffer, 13167,
                                   11, nmax);
    }

    for (size_t m = 0; m < 330; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
