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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecSDG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDS.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_sdg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_sdg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 7793, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 270 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 7793, 7199, 540, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 9, 6, 7,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 18, 6, 7,
                                                             ncols, fj, i * nprim_b + j, fq);

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

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 141, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 144, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 147, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 150, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 153, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 156, 0, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 159, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 162, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 165, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 168, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 171, 0, 6, 10, 11,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 180, 0, 6, 11, 12,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 189, 0, 6, 12, 13,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 198, 0, 6, 13, 14,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 207, 0, 6, 14, 15,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 216, 0, 6, 15, 16,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 225, 0, 6, 19, 20,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 234, 0, 6, 20, 21,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 243, 0, 6, 21, 22,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 252, 0, 6, 22, 23,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 261, 0, 6, 23, 24,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 270, 0, 6, 24, 25,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 279, 0, 3, 6, 27,
                                                                       30, 69, 75, 171, 180,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 297, 0, 3, 6, 30,
                                                                       33, 75, 81, 180, 189,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 315, 0, 3, 6, 33,
                                                                       36, 81, 87, 189, 198,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 333, 0, 3, 6, 36,
                                                                       39, 87, 93, 198, 207,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 351, 0, 3, 6, 39,
                                                                       42, 93, 99, 207, 216,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 369, 0, 3, 6, 48,
                                                                       51, 105, 111, 225, 234,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 387, 0, 3, 6, 51,
                                                                       54, 111, 117, 234, 243,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 405, 0, 3, 6, 54,
                                                                       57, 117, 123, 243, 252,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 423, 0, 3, 6, 57,
                                                                       60, 123, 129, 252, 261,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 441, 0, 3, 6, 60,
                                                                       63, 129, 135, 261, 270,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 459, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 462, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 465, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 468, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 471, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 474, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 477, 6, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 480, 6, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 483, 6, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 486, 6, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 489, 6, 25, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 492, 6, 26, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 495, 6, 12, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 504, 6, 13, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 513, 6, 14, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 522, 6, 15, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 531, 6, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 540, 6, 21, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 549, 6, 22, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 558, 6, 23, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 567, 6, 24, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 576, 6, 25, 66,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 585, 6, 33, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 603, 6, 36, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 621, 6, 39, 93,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 639, 6, 42, 99,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 657, 6, 54, 117,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 675, 6, 57, 123,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 693, 6, 60, 129,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 711, 6, 63, 135,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 729, 6, 12, 141,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 738, 6, 13, 144,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 747, 6, 14, 147,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 756, 6, 15, 150,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 765, 6, 16, 153,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 774, 6, 21, 156,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 783, 6, 22, 159,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 792, 6, 23, 162,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 801, 6, 24, 165,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 810, 6, 25, 168,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 819, 6, 33, 141,
                                                                       189, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 846, 6, 36, 144,
                                                                       198, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 873, 6, 39, 147,
                                                                       207, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 900, 6, 42, 150,
                                                                       216, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 927, 6, 54, 156,
                                                                       243, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 954, 6, 57, 159,
                                                                       252, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 981, 6, 60, 162,
                                                                       261, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1008, 6, 63, 165,
                                                                       270, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1035, 3, 6, 81,
                                                                       819, 189, 846, 315, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1089, 3, 6, 87,
                                                                       846, 198, 873, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1143, 3, 6, 93,
                                                                       873, 207, 900, 351, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1197, 3, 6, 117,
                                                                       927, 243, 954, 405, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1251, 3, 6, 123,
                                                                       954, 252, 981, 423, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1305, 3, 6, 129,
                                                                       981, 261, 1008, 441,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1359, 6, 10, 11,
                                                                       459, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1365, 6, 11, 12,
                                                                       462, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1371, 6, 12, 13,
                                                                       465, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1377, 6, 13, 14,
                                                                       468, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1383, 6, 14, 15,
                                                                       471, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1389, 6, 15, 16,
                                                                       474, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1395, 6, 19, 20,
                                                                       477, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1401, 6, 20, 21,
                                                                       480, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1407, 6, 21, 22,
                                                                       483, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1413, 6, 22, 23,
                                                                       486, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1419, 6, 23, 24,
                                                                       489, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1425, 6, 24, 25,
                                                                       492, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1431, 3, 6, 1359,
                                                                       459, 1365, 495, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1449, 3, 6, 1365,
                                                                       462, 1371, 504, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1467, 3, 6, 1371,
                                                                       465, 1377, 513, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1485, 3, 6, 1377,
                                                                       468, 1383, 522, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1503, 3, 6, 1383,
                                                                       471, 1389, 531, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1521, 3, 6, 1395,
                                                                       477, 1401, 540, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1539, 3, 6, 1401,
                                                                       480, 1407, 549, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1557, 3, 6, 1407,
                                                                       483, 1413, 558, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1575, 3, 6, 1413,
                                                                       486, 1419, 567, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1593, 3, 6, 1419,
                                                                       489, 1425, 576, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1611, 3, 6, 1431,
                                                                       495, 1449, 69, 75, 585,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1647, 3, 6, 1449,
                                                                       504, 1467, 75, 81, 603,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1683, 3, 6, 1467,
                                                                       513, 1485, 81, 87, 621,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1719, 3, 6, 1485,
                                                                       522, 1503, 87, 93, 639,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1755, 3, 6, 1521,
                                                                       540, 1539, 105, 111, 657,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1791, 3, 6, 1539,
                                                                       549, 1557, 111, 117, 675,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1827, 3, 6, 1557,
                                                                       558, 1575, 117, 123, 693,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 1863, 3, 6, 1575,
                                                                       567, 1593, 123, 129, 711,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1899, 0, 6, 1359,
                                                                       459, 1365, 729, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1917, 0, 6, 1365,
                                                                       462, 1371, 738, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1935, 0, 6, 1371,
                                                                       465, 1377, 747, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1953, 0, 6, 1377,
                                                                       468, 1383, 756, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1971, 0, 6, 1383,
                                                                       471, 1389, 765, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1989, 0, 6, 1395,
                                                                       477, 1401, 774, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2007, 0, 6, 1401,
                                                                       480, 1407, 783, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2025, 0, 6, 1407,
                                                                       483, 1413, 792, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2043, 0, 6, 1413,
                                                                       486, 1419, 801, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2061, 0, 6, 1419,
                                                                       489, 1425, 810, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2079, 0, 3, 6,
                                                                       1431, 495, 1449, 1899,
                                                                       729, 1917, 171, 180, 819,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2133, 0, 3, 6,
                                                                       1449, 504, 1467, 1917,
                                                                       738, 1935, 180, 189, 846,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2187, 0, 3, 6,
                                                                       1467, 513, 1485, 1935,
                                                                       747, 1953, 189, 198, 873,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2241, 0, 3, 6,
                                                                       1485, 522, 1503, 1953,
                                                                       756, 1971, 198, 207, 900,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2295, 0, 3, 6,
                                                                       1521, 540, 1539, 1989,
                                                                       774, 2007, 225, 234, 927,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2349, 0, 3, 6,
                                                                       1539, 549, 1557, 2007,
                                                                       783, 2025, 234, 243, 954,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2403, 0, 3, 6,
                                                                       1557, 558, 1575, 2025,
                                                                       792, 2043, 243, 252, 981,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2457, 0, 3, 6,
                                                                       1575, 567, 1593, 2043,
                                                                       801, 2061, 252, 261, 1008,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 2511, 0, 3, 6,
                                                                       1611, 585, 1647, 2079,
                                                                       819, 2133, 279, 297, 1035,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 2619, 0, 3, 6,
                                                                       1647, 603, 1683, 2133,
                                                                       846, 2187, 297, 315, 1089,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 2727, 0, 3, 6,
                                                                       1683, 621, 1719, 2187,
                                                                       873, 2241, 315, 333, 1143,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 2835, 0, 3, 6,
                                                                       1755, 657, 1791, 2295,
                                                                       927, 2349, 369, 387, 1197,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 2943, 0, 3, 6,
                                                                       1791, 675, 1827, 2349,
                                                                       954, 2403, 387, 405, 1251,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3051, 0, 3, 6,
                                                                       1827, 693, 1863, 2403,
                                                                       981, 2457, 405, 423, 1305,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3159, 6, 459, 462,
                                                                       1371, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3169, 6, 462, 465,
                                                                       1377, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3179, 6, 465, 468,
                                                                       1383, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3189, 6, 468, 471,
                                                                       1389, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3199, 6, 477, 480,
                                                                       1407, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3209, 6, 480, 483,
                                                                       1413, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3219, 6, 483, 486,
                                                                       1419, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3229, 6, 486, 489,
                                                                       1425, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3239, 3, 6, 3159,
                                                                       1371, 3169, 1467, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3269, 3, 6, 3169,
                                                                       1377, 3179, 1485, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3299, 3, 6, 3179,
                                                                       1383, 3189, 1503, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3329, 3, 6, 3199,
                                                                       1407, 3209, 1557, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3359, 3, 6, 3209,
                                                                       1413, 3219, 1575, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3389, 3, 6, 3219,
                                                                       1419, 3229, 1593, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3419, 3, 6, 3239,
                                                                       1467, 3269, 585, 603,
                                                                       1683, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3479, 3, 6, 3269,
                                                                       1485, 3299, 603, 621,
                                                                       1719, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3539, 3, 6, 3329,
                                                                       1557, 3359, 657, 675,
                                                                       1827, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 3599, 3, 6, 3359,
                                                                       1575, 3389, 675, 693,
                                                                       1863, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3659, 0, 6, 3159,
                                                                       1371, 3169, 729, 738,
                                                                       1935, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3689, 0, 6, 3169,
                                                                       1377, 3179, 738, 747,
                                                                       1953, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3719, 0, 6, 3179,
                                                                       1383, 3189, 747, 756,
                                                                       1971, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3749, 0, 6, 3199,
                                                                       1407, 3209, 774, 783,
                                                                       2025, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3779, 0, 6, 3209,
                                                                       1413, 3219, 783, 792,
                                                                       2043, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3809, 0, 6, 3219,
                                                                       1419, 3229, 792, 801,
                                                                       2061, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 3839, 0, 3, 6,
                                                                       3239, 1467, 3269, 3659,
                                                                       1935, 3689, 819, 846,
                                                                       2187, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 3929, 0, 3, 6,
                                                                       3269, 1485, 3299, 3689,
                                                                       1953, 3719, 846, 873,
                                                                       2241, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4019, 0, 3, 6,
                                                                       3329, 1557, 3359, 3749,
                                                                       2025, 3779, 927, 954,
                                                                       2403, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4109, 0, 3, 6,
                                                                       3359, 1575, 3389, 3779,
                                                                       2043, 3809, 954, 981,
                                                                       2457, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 4199, 0, 3, 6,
                                                                       3419, 1683, 3479, 3839,
                                                                       2187, 3929, 1035, 1089,
                                                                       2727, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 4379, 0, 3, 6,
                                                                       3539, 1827, 3599, 4019,
                                                                       2403, 4109, 1197, 1251,
                                                                       3051, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4559, 6, 1359,
                                                                       1365, 3159, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4574, 6, 1365,
                                                                       1371, 3169, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4589, 6, 1371,
                                                                       1377, 3179, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4604, 6, 1377,
                                                                       1383, 3189, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4619, 6, 1395,
                                                                       1401, 3199, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4634, 6, 1401,
                                                                       1407, 3209, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4649, 6, 1407,
                                                                       1413, 3219, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 4664, 6, 1413,
                                                                       1419, 3229, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4679, 3, 6, 4559,
                                                                       3159, 4574, 1431, 1449,
                                                                       3239, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4724, 3, 6, 4574,
                                                                       3169, 4589, 1449, 1467,
                                                                       3269, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4769, 3, 6, 4589,
                                                                       3179, 4604, 1467, 1485,
                                                                       3299, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4814, 3, 6, 4619,
                                                                       3199, 4634, 1521, 1539,
                                                                       3329, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4859, 3, 6, 4634,
                                                                       3209, 4649, 1539, 1557,
                                                                       3359, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 4904, 3, 6, 4649,
                                                                       3219, 4664, 1557, 1575,
                                                                       3389, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 4949, 3, 6, 4679,
                                                                       3239, 4724, 1611, 1647,
                                                                       3419, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5039, 3, 6, 4724,
                                                                       3269, 4769, 1647, 1683,
                                                                       3479, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5129, 3, 6, 4814,
                                                                       3329, 4859, 1755, 1791,
                                                                       3539, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 5219, 3, 6, 4859,
                                                                       3359, 4904, 1791, 1827,
                                                                       3599, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5309, 0, 6, 4559,
                                                                       3159, 4574, 1899, 1917,
                                                                       3659, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5354, 0, 6, 4574,
                                                                       3169, 4589, 1917, 1935,
                                                                       3689, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5399, 0, 6, 4589,
                                                                       3179, 4604, 1935, 1953,
                                                                       3719, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5444, 0, 6, 4619,
                                                                       3199, 4634, 1989, 2007,
                                                                       3749, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5489, 0, 6, 4634,
                                                                       3209, 4649, 2007, 2025,
                                                                       3779, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 5534, 0, 6, 4649,
                                                                       3219, 4664, 2025, 2043,
                                                                       3809, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 5579, 0, 3, 6,
                                                                       4679, 3239, 4724, 5309,
                                                                       3659, 5354, 2079, 2133,
                                                                       3839, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 5714, 0, 3, 6,
                                                                       4724, 3269, 4769, 5354,
                                                                       3689, 5399, 2133, 2187,
                                                                       3929, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 5849, 0, 3, 6,
                                                                       4814, 3329, 4859, 5444,
                                                                       3749, 5489, 2295, 2349,
                                                                       4019, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 5984, 0, 3, 6,
                                                                       4859, 3359, 4904, 5489,
                                                                       3779, 5534, 2349, 2403,
                                                                       4109, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 6119, 0, 3, 6,
                                                                       4949, 3419, 5039, 5579,
                                                                       3839, 5714, 2511, 2619,
                                                                       4199, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 6389, 0, 3, 6,
                                                                       5129, 3539, 5219, 5849,
                                                                       4019, 5984, 2835, 2943,
                                                                       4379, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_s_x(buffer, 6659, 6389, 1, 90, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 6749, 6389, 1, 90, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 6839, 6389, 1, 90, ncols, alpha);

                    simdgeo::geom_s_x(buffer, 6929, 6119, 1, 90, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 7019, 6119, 1, 90, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 7109, 6119, 1, 90, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 7199, 6659, 540, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 7739, 7199, 6, 1, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 7739, 9, nmax);

        simdtrf::transform_g_inner(buffer, 7739, 7289, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 45 * nvalues + n * npairs, nvalues, buffer, 7739, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 7739, 7379, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 90 * nvalues + n * npairs, nvalues, buffer, 7739, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 7739, 7469, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 135 * nvalues + n * npairs, nvalues, buffer, 7739, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 7739, 7559, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 180 * nvalues + n * npairs, nvalues, buffer, 7739, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 7739, 7649, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 225 * nvalues + n * npairs, nvalues, buffer, 7739, 9,
                                   nmax);
    }

    for (size_t m = 0; m < 270; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
