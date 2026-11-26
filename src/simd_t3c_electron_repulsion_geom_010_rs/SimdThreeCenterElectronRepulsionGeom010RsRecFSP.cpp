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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecFSP.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_fsp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_fsp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 1279, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 126 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 1279, 1069, 180, dimensions);

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
                                                            4, 5}, ncols, fj, i * nprim_b + j,
                                                            fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 15, 6, {1, 2, 3, 4,
                                                        5}, ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 3, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 60, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 63, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 66, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 69, 0, 6, 10, 11,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 78, 0, 6, 11, 12,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 87, 0, 6, 12, 13,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 96, 0, 6, 16, 17,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 105, 0, 6, 17, 18,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 114, 0, 6, 18, 19,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 123, 0, 6, 10, 11,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 129, 0, 6, 11, 12,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 135, 0, 6, 12, 13,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 141, 0, 6, 16, 17,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 147, 0, 6, 17, 18,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 153, 0, 6, 18, 19,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 159, 0, 3, 6, 45,
                                                                       48, 69, 78, 123, 129,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 177, 0, 3, 6, 48,
                                                                       51, 78, 87, 129, 135,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 195, 0, 3, 6, 57,
                                                                       60, 96, 105, 141, 147,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 213, 0, 3, 6, 60,
                                                                       63, 105, 114, 147, 153,
                                                                       ncols, gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 231, 0, 6, 45, 48,
                                                                       123, 129, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 241, 0, 6, 48, 51,
                                                                       129, 135, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 251, 0, 6, 57, 60,
                                                                       141, 147, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 261, 0, 6, 60, 63,
                                                                       147, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 271, 0, 3, 6, 123,
                                                                       129, 159, 177, 231, 241,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 301, 0, 3, 6, 141,
                                                                       147, 195, 213, 251, 261,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 331, 6, 21, 45,
                                                                       69, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 358, 6, 24, 48,
                                                                       78, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 385, 6, 27, 51,
                                                                       87, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 412, 6, 33, 57,
                                                                       96, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 439, 6, 36, 60,
                                                                       105, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 466, 6, 39, 63,
                                                                       114, ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 493, 0, 6, 331,
                                                                       69, 358, 123, 159, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 547, 0, 6, 358,
                                                                       78, 385, 129, 177, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 601, 0, 6, 412,
                                                                       96, 439, 141, 195, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 655, 0, 6, 439,
                                                                       105, 466, 147, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 709, 0, 6, 493,
                                                                       159, 547, 231, 271, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 799, 0, 6, 601,
                                                                       195, 655, 251, 301, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 889, 799, 10, 3, ncols, beta);

                    simdgeo::geom_s_y(buffer, 919, 799, 10, 3, ncols, beta);

                    simdgeo::geom_s_z(buffer, 949, 799, 10, 3, ncols, beta);

                    simdgeo::geom_s_x(buffer, 979, 709, 10, 3, ncols, beta);

                    simdgeo::geom_s_y(buffer, 1009, 709, 10, 3, ncols, beta);

                    simdgeo::geom_s_z(buffer, 1039, 709, 10, 3, ncols, beta);

                    simdfunc::contract_primitives(buffer, 1069, 889, 180, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 1249, 1069, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 1249, 3, nmax);

        simdtrf::transform_p_inner(buffer, 1249, 1099, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 21 * nvalues + n * npairs, nvalues, buffer, 1249, 3,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 1249, 1129, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 42 * nvalues + n * npairs, nvalues, buffer, 1249, 3,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 1249, 1159, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 1249, 3,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 1249, 1189, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 84 * nvalues + n * npairs, nvalues, buffer, 1249, 3,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 1249, 1219, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 105 * nvalues + n * npairs, nvalues, buffer, 1249, 3,
                                   nmax);
    }

    for (size_t m = 0; m < 126; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
