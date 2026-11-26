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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecSFP.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_sfp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_sfp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 1273, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 1273, 1063, 180, dimensions);

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

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 45, 3, 6, 10, 11,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 51, 3, 6, 11, 12,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 57, 3, 6, 12, 13,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 63, 3, 6, 16, 17,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 69, 3, 6, 17, 18,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 75, 3, 6, 18, 19,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 81, 3, 6, 21, 24,
                                                                       45, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 91, 3, 6, 24, 27,
                                                                       51, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 101, 3, 6, 33, 36,
                                                                       63, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 111, 3, 6, 36, 39,
                                                                       69, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 121, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 124, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 127, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 130, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 133, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 136, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 139, 0, 6, 10, 11,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 148, 0, 6, 11, 12,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 157, 0, 6, 12, 13,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 166, 0, 6, 16, 17,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 175, 0, 6, 17, 18,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 184, 0, 6, 18, 19,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 193, 0, 3, 6, 21,
                                                                       24, 45, 51, 139, 148,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 211, 0, 3, 6, 24,
                                                                       27, 51, 57, 148, 157,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 229, 0, 3, 6, 33,
                                                                       36, 63, 69, 166, 175,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 247, 0, 3, 6, 36,
                                                                       39, 69, 75, 175, 184,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 265, 0, 3, 6, 45,
                                                                       51, 81, 91, 193, 211,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 295, 0, 3, 6, 63,
                                                                       69, 101, 111, 229, 247,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 325, 6, 21, 121,
                                                                       139, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 352, 6, 24, 124,
                                                                       148, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 379, 6, 27, 127,
                                                                       157, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 406, 6, 33, 130,
                                                                       166, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 433, 6, 36, 133,
                                                                       175, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 460, 6, 39, 136,
                                                                       184, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 487, 3, 6, 45,
                                                                       325, 139, 352, 193, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 541, 3, 6, 51,
                                                                       352, 148, 379, 211, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 595, 3, 6, 63,
                                                                       406, 166, 433, 229, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 649, 3, 6, 69,
                                                                       433, 175, 460, 247, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 703, 3, 6, 81,
                                                                       487, 193, 541, 265, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 793, 3, 6, 101,
                                                                       595, 229, 649, 295, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 883, 793, 1, 30, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 913, 793, 1, 30, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 943, 793, 1, 30, ncols, alpha);

                    simdgeo::geom_s_x(buffer, 973, 703, 1, 30, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 1003, 703, 1, 30, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 1033, 703, 1, 30, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 1063, 883, 180, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 1243, 1063, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 1243, 3, nmax);

        simdtrf::transform_p_inner(buffer, 1243, 1093, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 21 * nvalues + n * npairs, nvalues, buffer, 1243, 3,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 1243, 1123, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 42 * nvalues + n * npairs, nvalues, buffer, 1243, 3,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 1243, 1153, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 1243, 3,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 1243, 1183, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 84 * nvalues + n * npairs, nvalues, buffer, 1243, 3,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 1243, 1213, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 105 * nvalues + n * npairs, nvalues, buffer, 1243, 3,
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
