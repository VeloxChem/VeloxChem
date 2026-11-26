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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecPPP.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdGeometryD1.hpp"
#include "SimdGeometryP1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferGeom010XPP.hpp"
#include "SimdTransferGeom010YPP.hpp"
#include "SimdTransferGeom010ZPP.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_ppp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_ppp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 909, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 162 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 909, 360, 342, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto beta = b_exps[j];

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pb(buffer, coordinates, 0, nmax, fb);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 3, n, nmax, fc);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_t3c_erf_boys_function(buffer, coordinates, 6, 3, {1, 2, 3,
                                                            4}, ncols, fj, i * nprim_b + j, fq,
                                                            omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 11, 3, {1, 2, 3, 4},
                                                        ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 16, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 19, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 34, 0, 3, 7, 8,
                                                                       16, 19, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 40, 0, 3, 8, 9,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 46, 0, 3, 12, 13,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 52, 0, 3, 13, 14,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 58, 0, 3, 16, 19,
                                                                       34, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 68, 0, 3, 25, 28,
                                                                       46, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 78, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 81, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 84, 3, 7, 16,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 93, 3, 12, 25,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 102, 3, 16, 34,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 120, 3, 25, 46,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 138, 3, 34, 58,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 168, 3, 46, 68,
                                                                       ncols, p, q);

                    simdgeo::geom_p_x(buffer, 198, 78, 102, 1, 3, ncols, beta);

                    simdgeo::geom_p_y(buffer, 207, 78, 102, 1, 3, ncols, beta);

                    simdgeo::geom_p_z(buffer, 216, 78, 102, 1, 3, ncols, beta);

                    simdgeo::geom_p_x(buffer, 225, 81, 120, 1, 3, ncols, beta);

                    simdgeo::geom_p_y(buffer, 234, 81, 120, 1, 3, ncols, beta);

                    simdgeo::geom_p_z(buffer, 243, 81, 120, 1, 3, ncols, beta);

                    simdgeo::geom_d_x(buffer, 252, 84, 138, 1, 3, ncols, beta);

                    simdgeo::geom_d_y(buffer, 270, 84, 138, 1, 3, ncols, beta);

                    simdgeo::geom_d_z(buffer, 288, 84, 138, 1, 3, ncols, beta);

                    simdgeo::geom_d_x(buffer, 306, 93, 168, 1, 3, ncols, beta);

                    simdgeo::geom_d_y(buffer, 324, 93, 168, 1, 3, ncols, beta);

                    simdgeo::geom_d_z(buffer, 342, 93, 168, 1, 3, ncols, beta);

                    simdfunc::contract_primitives(buffer, 360, 198, 9, ncols);

                    simdfunc::contract_primitives(buffer, 378, 207, 9, ncols);

                    simdfunc::contract_primitives(buffer, 396, 216, 9, ncols);

                    simdfunc::contract_primitives(buffer, 414, 84, 9, ncols);

                    simdfunc::contract_primitives(buffer, 432, 225, 9, ncols);

                    simdfunc::contract_primitives(buffer, 450, 234, 9, ncols);

                    simdfunc::contract_primitives(buffer, 468, 243, 9, ncols);

                    simdfunc::contract_primitives(buffer, 486, 93, 9, ncols);

                    simdfunc::contract_primitives(buffer, 504, 252, 18, ncols);

                    simdfunc::contract_primitives(buffer, 540, 270, 18, ncols);

                    simdfunc::contract_primitives(buffer, 576, 288, 18, ncols);

                    simdfunc::contract_primitives(buffer, 612, 306, 18, ncols);

                    simdfunc::contract_primitives(buffer, 648, 324, 18, ncols);

                    simdfunc::contract_primitives(buffer, 684, 342, 18, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 369, 360, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 387, 378, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 405, 396, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 423, 414, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 441, 432, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 459, 450, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 477, 468, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 495, 486, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 522, 504, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 558, 540, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 594, 576, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 630, 612, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 666, 648, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 702, 684, 6, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 720, 369, 423, 522, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 747, 387, 423, 558, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 774, 405, 423, 594, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 801, 441, 495, 630, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 828, 459, 495, 666, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 855, 477, 495, 702, 3, nmax);

        simdtrf::transform_p_inner(buffer, 882, 801, 3, 3, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 882, 9, nmax);

        simdtrf::transform_p_inner(buffer, 882, 828, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 27 * nvalues + n * npairs, nvalues, buffer, 882, 9,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 882, 855, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 54 * nvalues + n * npairs, nvalues, buffer, 882, 9,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 882, 720, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 81 * nvalues + n * npairs, nvalues, buffer, 882, 9,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 882, 747, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 108 * nvalues + n * npairs, nvalues, buffer, 882, 9,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 882, 774, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 135 * nvalues + n * npairs, nvalues, buffer, 882, 9,
                                   nmax);
    }

    for (size_t m = 0; m < 162; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
