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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecSDH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDS.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
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
compute_rs_geom_100_sdh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_sdh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 69, 3, 6, 10, 11,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 75, 3, 6, 11, 12,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 81, 3, 6, 12, 13,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 87, 3, 6, 13, 14,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 93, 3, 6, 14, 15,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 99, 3, 6, 15, 16,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 105, 3, 6, 19, 20,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 111, 3, 6, 20, 21,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 117, 3, 6, 21, 22,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 123, 3, 6, 22, 23,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 129, 3, 6, 23, 24,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 135, 3, 6, 24, 25,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 141, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 144, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 147, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 150, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 153, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 156, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 159, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 162, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 165, 0, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 168, 0, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 171, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 174, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 177, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 180, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 183, 0, 6, 10, 11,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 192, 0, 6, 11, 12,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 201, 0, 6, 12, 13,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 210, 0, 6, 13, 14,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 219, 0, 6, 14, 15,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 228, 0, 6, 15, 16,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 237, 0, 6, 19, 20,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 246, 0, 6, 20, 21,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 255, 0, 6, 21, 22,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 264, 0, 6, 22, 23,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 273, 0, 6, 23, 24,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 282, 0, 6, 24, 25,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 291, 0, 3, 6, 27,
                                                                       30, 69, 75, 183, 192,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 309, 0, 3, 6, 30,
                                                                       33, 75, 81, 192, 201,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 327, 0, 3, 6, 33,
                                                                       36, 81, 87, 201, 210,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 345, 0, 3, 6, 36,
                                                                       39, 87, 93, 210, 219,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 363, 0, 3, 6, 39,
                                                                       42, 93, 99, 219, 228,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 381, 0, 3, 6, 48,
                                                                       51, 105, 111, 237, 246,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 399, 0, 3, 6, 51,
                                                                       54, 111, 117, 246, 255,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 417, 0, 3, 6, 54,
                                                                       57, 117, 123, 255, 264,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 435, 0, 3, 6, 57,
                                                                       60, 123, 129, 264, 273,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 453, 0, 3, 6, 60,
                                                                       63, 129, 135, 273, 282,
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

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 609, 6, 27, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 627, 6, 30, 75,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 645, 6, 33, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 663, 6, 36, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 681, 6, 39, 93,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 699, 6, 42, 99,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 717, 6, 48, 105,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 735, 6, 51, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 753, 6, 54, 117,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 771, 6, 57, 123,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 789, 6, 60, 129,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 807, 6, 63, 135,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 825, 6, 10, 141,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 834, 6, 11, 144,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 843, 6, 12, 147,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 852, 6, 13, 150,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 861, 6, 14, 153,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 870, 6, 15, 156,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 879, 6, 16, 159,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 888, 6, 19, 162,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 897, 6, 20, 165,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 906, 6, 21, 168,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 915, 6, 22, 171,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 924, 6, 23, 174,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 933, 6, 24, 177,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 942, 6, 25, 180,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 951, 6, 27, 141,
                                                                       183, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 978, 6, 30, 144,
                                                                       192, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1005, 6, 33, 147,
                                                                       201, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1032, 6, 36, 150,
                                                                       210, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1059, 6, 39, 153,
                                                                       219, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1086, 6, 42, 156,
                                                                       228, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1113, 6, 48, 162,
                                                                       237, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1140, 6, 51, 165,
                                                                       246, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1167, 6, 54, 168,
                                                                       255, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1194, 6, 57, 171,
                                                                       264, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1221, 6, 60, 174,
                                                                       273, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1248, 6, 63, 177,
                                                                       282, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1275, 3, 6, 69,
                                                                       951, 183, 978, 291, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1329, 3, 6, 75,
                                                                       978, 192, 1005, 309,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1383, 3, 6, 81,
                                                                       1005, 201, 1032, 327,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1437, 3, 6, 87,
                                                                       1032, 210, 1059, 345,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1491, 3, 6, 93,
                                                                       1059, 219, 1086, 363,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1545, 3, 6, 105,
                                                                       1113, 237, 1140, 381,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1599, 3, 6, 111,
                                                                       1140, 246, 1167, 399,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1653, 3, 6, 117,
                                                                       1167, 255, 1194, 417,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1707, 3, 6, 123,
                                                                       1194, 264, 1221, 435,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1761, 3, 6, 129,
                                                                       1221, 273, 1248, 453,
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

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2067, 3, 6, 1887,
                                                                       519, 1905, 69, 75, 645,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2103, 3, 6, 1905,
                                                                       528, 1923, 75, 81, 663,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2139, 3, 6, 1923,
                                                                       537, 1941, 81, 87, 681,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2175, 3, 6, 1941,
                                                                       546, 1959, 87, 93, 699,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2211, 3, 6, 1977,
                                                                       564, 1995, 105, 111, 753,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2247, 3, 6, 1995,
                                                                       573, 2013, 111, 117, 771,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2283, 3, 6, 2013,
                                                                       582, 2031, 117, 123, 789,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2319, 3, 6, 2031,
                                                                       591, 2049, 123, 129, 807,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2355, 0, 6, 1815,
                                                                       477, 1821, 843, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2373, 0, 6, 1821,
                                                                       480, 1827, 852, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2391, 0, 6, 1827,
                                                                       483, 1833, 861, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2409, 0, 6, 1833,
                                                                       486, 1839, 870, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2427, 0, 6, 1839,
                                                                       489, 1845, 879, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2445, 0, 6, 1851,
                                                                       501, 1857, 906, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2463, 0, 6, 1857,
                                                                       504, 1863, 915, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2481, 0, 6, 1863,
                                                                       507, 1869, 924, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2499, 0, 6, 1869,
                                                                       510, 1875, 933, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2517, 0, 6, 1875,
                                                                       513, 1881, 942, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2535, 0, 3, 6,
                                                                       1887, 519, 1905, 2355,
                                                                       843, 2373, 183, 192, 1005,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2589, 0, 3, 6,
                                                                       1905, 528, 1923, 2373,
                                                                       852, 2391, 192, 201, 1032,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2643, 0, 3, 6,
                                                                       1923, 537, 1941, 2391,
                                                                       861, 2409, 201, 210, 1059,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2697, 0, 3, 6,
                                                                       1941, 546, 1959, 2409,
                                                                       870, 2427, 210, 219, 1086,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2751, 0, 3, 6,
                                                                       1977, 564, 1995, 2445,
                                                                       906, 2463, 237, 246, 1167,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2805, 0, 3, 6,
                                                                       1995, 573, 2013, 2463,
                                                                       915, 2481, 246, 255, 1194,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2859, 0, 3, 6,
                                                                       2013, 582, 2031, 2481,
                                                                       924, 2499, 255, 264, 1221,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2913, 0, 3, 6,
                                                                       2031, 591, 2049, 2499,
                                                                       933, 2517, 264, 273, 1248,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 2967, 0, 3, 6,
                                                                       2067, 645, 2103, 2535,
                                                                       1005, 2589, 291, 309,
                                                                       1383, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3075, 0, 3, 6,
                                                                       2103, 663, 2139, 2589,
                                                                       1032, 2643, 309, 327,
                                                                       1437, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3183, 0, 3, 6,
                                                                       2139, 681, 2175, 2643,
                                                                       1059, 2697, 327, 345,
                                                                       1491, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3291, 0, 3, 6,
                                                                       2211, 753, 2247, 2751,
                                                                       1167, 2805, 381, 399,
                                                                       1653, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3399, 0, 3, 6,
                                                                       2247, 771, 2283, 2805,
                                                                       1194, 2859, 399, 417,
                                                                       1707, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3507, 0, 3, 6,
                                                                       2283, 789, 2319, 2859,
                                                                       1221, 2913, 417, 435,
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

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4035, 3, 6, 3735,
                                                                       1887, 3765, 609, 627,
                                                                       2067, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4095, 3, 6, 3765,
                                                                       1905, 3795, 627, 645,
                                                                       2103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4155, 3, 6, 3795,
                                                                       1923, 3825, 645, 663,
                                                                       2139, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4215, 3, 6, 3825,
                                                                       1941, 3855, 663, 681,
                                                                       2175, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4275, 3, 6, 3885,
                                                                       1977, 3915, 717, 735,
                                                                       2211, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4335, 3, 6, 3915,
                                                                       1995, 3945, 735, 753,
                                                                       2247, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4395, 3, 6, 3945,
                                                                       2013, 3975, 753, 771,
                                                                       2283, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4455, 3, 6, 3975,
                                                                       2031, 4005, 771, 789,
                                                                       2319, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4515, 0, 6, 3615,
                                                                       1815, 3625, 825, 834,
                                                                       2355, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4545, 0, 6, 3625,
                                                                       1821, 3635, 834, 843,
                                                                       2373, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4575, 0, 6, 3635,
                                                                       1827, 3645, 843, 852,
                                                                       2391, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4605, 0, 6, 3645,
                                                                       1833, 3655, 852, 861,
                                                                       2409, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4635, 0, 6, 3655,
                                                                       1839, 3665, 861, 870,
                                                                       2427, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4665, 0, 6, 3675,
                                                                       1851, 3685, 888, 897,
                                                                       2445, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4695, 0, 6, 3685,
                                                                       1857, 3695, 897, 906,
                                                                       2463, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4725, 0, 6, 3695,
                                                                       1863, 3705, 906, 915,
                                                                       2481, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4755, 0, 6, 3705,
                                                                       1869, 3715, 915, 924,
                                                                       2499, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4785, 0, 6, 3715,
                                                                       1875, 3725, 924, 933,
                                                                       2517, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4815, 0, 3, 6,
                                                                       3735, 1887, 3765, 4515,
                                                                       2355, 4545, 951, 978,
                                                                       2535, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4905, 0, 3, 6,
                                                                       3765, 1905, 3795, 4545,
                                                                       2373, 4575, 978, 1005,
                                                                       2589, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4995, 0, 3, 6,
                                                                       3795, 1923, 3825, 4575,
                                                                       2391, 4605, 1005, 1032,
                                                                       2643, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5085, 0, 3, 6,
                                                                       3825, 1941, 3855, 4605,
                                                                       2409, 4635, 1032, 1059,
                                                                       2697, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5175, 0, 3, 6,
                                                                       3885, 1977, 3915, 4665,
                                                                       2445, 4695, 1113, 1140,
                                                                       2751, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5265, 0, 3, 6,
                                                                       3915, 1995, 3945, 4695,
                                                                       2463, 4725, 1140, 1167,
                                                                       2805, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5355, 0, 3, 6,
                                                                       3945, 2013, 3975, 4725,
                                                                       2481, 4755, 1167, 1194,
                                                                       2859, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5445, 0, 3, 6,
                                                                       3975, 2031, 4005, 4755,
                                                                       2499, 4785, 1194, 1221,
                                                                       2913, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 5535, 0, 3, 6,
                                                                       4035, 2067, 4095, 4815,
                                                                       2535, 4905, 1275, 1329,
                                                                       2967, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 5715, 0, 3, 6,
                                                                       4095, 2103, 4155, 4905,
                                                                       2589, 4995, 1329, 1383,
                                                                       3075, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 5895, 0, 3, 6,
                                                                       4155, 2139, 4215, 4995,
                                                                       2643, 5085, 1383, 1437,
                                                                       3183, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 6075, 0, 3, 6,
                                                                       4275, 2211, 4335, 5175,
                                                                       2751, 5265, 1545, 1599,
                                                                       3291, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 6255, 0, 3, 6,
                                                                       4335, 2247, 4395, 5265,
                                                                       2805, 5355, 1599, 1653,
                                                                       3399, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 6435, 0, 3, 6,
                                                                       4395, 2283, 4455, 5355,
                                                                       2859, 5445, 1653, 1707,
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

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7005, 3, 6, 6735,
                                                                       3795, 6780, 2067, 2103,
                                                                       4155, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7095, 3, 6, 6780,
                                                                       3825, 6825, 2103, 2139,
                                                                       4215, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7185, 3, 6, 6870,
                                                                       3945, 6915, 2211, 2247,
                                                                       4395, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7275, 3, 6, 6915,
                                                                       3975, 6960, 2247, 2283,
                                                                       4455, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7365, 0, 6, 6615,
                                                                       3635, 6630, 2355, 2373,
                                                                       4575, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7410, 0, 6, 6630,
                                                                       3645, 6645, 2373, 2391,
                                                                       4605, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7455, 0, 6, 6645,
                                                                       3655, 6660, 2391, 2409,
                                                                       4635, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7500, 0, 6, 6675,
                                                                       3695, 6690, 2445, 2463,
                                                                       4725, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7545, 0, 6, 6690,
                                                                       3705, 6705, 2463, 2481,
                                                                       4755, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7590, 0, 6, 6705,
                                                                       3715, 6720, 2481, 2499,
                                                                       4785, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 7635, 0, 3, 6,
                                                                       6735, 3795, 6780, 7365,
                                                                       4575, 7410, 2535, 2589,
                                                                       4995, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 7770, 0, 3, 6,
                                                                       6780, 3825, 6825, 7410,
                                                                       4605, 7455, 2589, 2643,
                                                                       5085, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 7905, 0, 3, 6,
                                                                       6870, 3945, 6915, 7500,
                                                                       4725, 7545, 2751, 2805,
                                                                       5355, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 8040, 0, 3, 6,
                                                                       6915, 3975, 6960, 7545,
                                                                       4755, 7590, 2805, 2859,
                                                                       5445, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 8175, 0, 3, 6,
                                                                       7005, 4155, 7095, 7635,
                                                                       4995, 7770, 2967, 3075,
                                                                       5895, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 8445, 0, 3, 6,
                                                                       7185, 4395, 7275, 7905,
                                                                       5355, 8040, 3291, 3399,
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

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 9261, 3, 6, 8883,
                                                                       6735, 8946, 4035, 4095,
                                                                       7005, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 9387, 3, 6, 8946,
                                                                       6780, 9009, 4095, 4155,
                                                                       7095, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 9513, 3, 6, 9072,
                                                                       6870, 9135, 4275, 4335,
                                                                       7185, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 9639, 3, 6, 9135,
                                                                       6915, 9198, 4335, 4395,
                                                                       7275, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9765, 0, 6, 8715,
                                                                       6615, 8736, 4515, 4545,
                                                                       7365, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9828, 0, 6, 8736,
                                                                       6630, 8757, 4545, 4575,
                                                                       7410, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9891, 0, 6, 8757,
                                                                       6645, 8778, 4575, 4605,
                                                                       7455, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 9954, 0, 6, 8799,
                                                                       6675, 8820, 4665, 4695,
                                                                       7500, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 10017, 0, 6, 8820,
                                                                       6690, 8841, 4695, 4725,
                                                                       7545, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 10080, 0, 6, 8841,
                                                                       6705, 8862, 4725, 4755,
                                                                       7590, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 10143, 0, 3, 6,
                                                                       8883, 6735, 8946, 9765,
                                                                       7365, 9828, 4815, 4905,
                                                                       7635, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 10332, 0, 3, 6,
                                                                       8946, 6780, 9009, 9828,
                                                                       7410, 9891, 4905, 4995,
                                                                       7770, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 10521, 0, 3, 6,
                                                                       9072, 6870, 9135, 9954,
                                                                       7500, 10017, 5175, 5265,
                                                                       7905, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 10710, 0, 3, 6,
                                                                       9135, 6915, 9198, 10017,
                                                                       7545, 10080, 5265, 5355,
                                                                       8040, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 10899, 0, 3, 6,
                                                                       9261, 7005, 9387, 10143,
                                                                       7635, 10332, 5535, 5715,
                                                                       8175, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 11277, 0, 3, 6,
                                                                       9513, 7185, 9639, 10521,
                                                                       7905, 10710, 6075, 6255,
                                                                       8445, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_s_x(buffer, 11655, 11277, 1, 126, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 11781, 11277, 1, 126, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 11907, 11277, 1, 126, ncols, alpha);

                    simdgeo::geom_s_x(buffer, 12033, 10899, 1, 126, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 12159, 10899, 1, 126, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 12285, 10899, 1, 126, ncols, alpha);

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
