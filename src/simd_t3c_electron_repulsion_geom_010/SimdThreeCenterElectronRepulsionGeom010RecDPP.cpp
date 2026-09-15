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


#include "SimdThreeCenterElectronRepulsionGeom010RecDPP.hpp"

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
#include "SimdGeometryF1.hpp"
#include "SimdGeometryP1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferGeom010XDP.hpp"
#include "SimdTransferGeom010XPD.hpp"
#include "SimdTransferGeom010XPP.hpp"
#include "SimdTransferGeom010YDP.hpp"
#include "SimdTransferGeom010YPD.hpp"
#include "SimdTransferGeom010YPP.hpp"
#include "SimdTransferGeom010ZDP.hpp"
#include "SimdTransferGeom010ZPD.hpp"
#include "SimdTransferGeom010ZPP.hpp"
#include "SimdTransferPP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_dpp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_dpp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 1235, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 135 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 1235, 353, 366, dimensions);

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 6, 3, {1, 2, 3, 4,
                                                        5}, ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 12, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 15, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 18, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 24, 0, 3, 7, 8,
                                                                       12, 15, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 30, 0, 3, 8, 9,
                                                                       15, 18, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 36, 0, 3, 9, 10,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 42, 0, 3, 12, 15,
                                                                       24, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 52, 0, 3, 15, 18,
                                                                       30, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 62, 0, 3, 24, 30,
                                                                       42, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 77, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 80, 3, 7, 12,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 89, 3, 12, 24,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 107, 3, 24, 42,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 137, 3, 42, 62,
                                                                       ncols, p, q);

                    simdgeo::geom_p_x(buffer, 182, 77, 89, 1, 3, ncols, beta);

                    simdgeo::geom_p_y(buffer, 191, 77, 89, 1, 3, ncols, beta);

                    simdgeo::geom_p_z(buffer, 200, 77, 89, 1, 3, ncols, beta);

                    simdgeo::geom_d_x(buffer, 209, 80, 107, 1, 3, ncols, beta);

                    simdgeo::geom_d_y(buffer, 227, 80, 107, 1, 3, ncols, beta);

                    simdgeo::geom_d_z(buffer, 245, 80, 107, 1, 3, ncols, beta);

                    simdgeo::geom_f_x(buffer, 263, 89, 137, 1, 3, ncols, beta);

                    simdgeo::geom_f_y(buffer, 293, 89, 137, 1, 3, ncols, beta);

                    simdgeo::geom_f_z(buffer, 323, 89, 137, 1, 3, ncols, beta);

                    simdfunc::contract_primitives(buffer, 353, 182, 9, ncols);

                    simdfunc::contract_primitives(buffer, 371, 191, 9, ncols);

                    simdfunc::contract_primitives(buffer, 389, 200, 9, ncols);

                    simdfunc::contract_primitives(buffer, 407, 80, 9, ncols);

                    simdfunc::contract_primitives(buffer, 425, 209, 18, ncols);

                    simdfunc::contract_primitives(buffer, 461, 227, 18, ncols);

                    simdfunc::contract_primitives(buffer, 497, 245, 18, ncols);

                    simdfunc::contract_primitives(buffer, 533, 89, 18, ncols);

                    simdfunc::contract_primitives(buffer, 569, 263, 30, ncols);

                    simdfunc::contract_primitives(buffer, 629, 293, 30, ncols);

                    simdfunc::contract_primitives(buffer, 689, 323, 30, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 362, 353, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 380, 371, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 398, 389, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 416, 407, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 443, 425, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 479, 461, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 515, 497, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 551, 533, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 599, 569, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 659, 629, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 719, 689, 10, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 749, 362, 416, 443, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 776, 380, 416, 479, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 803, 398, 416, 515, 3, nmax);

        simdtrf::compute_hrr_pp(buffer, coordinates, 830, 416, 551, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 857, 443, 551, 599, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 911, 479, 551, 659, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 965, 515, 551, 719, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dp_out_of_second(buffer, coordinates, 1019, 749, 830, 857,
                                                        3, nmax);

        simdtrf::compute_hrr_geom_010y_dp_out_of_second(buffer, coordinates, 1073, 776, 830, 911,
                                                        3, nmax);

        simdtrf::compute_hrr_geom_010z_dp_out_of_second(buffer, coordinates, 1127, 803, 830, 965,
                                                        3, nmax);

        simdtrf::transform_p_inner(buffer, 1181, 1019, 6, 3, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 1181, 9, nmax);

        simdtrf::transform_p_inner(buffer, 1181, 1073, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 45 * nvalues + n * npairs, nvalues, buffer, 1181, 9,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 1181, 1127, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 90 * nvalues + n * npairs, nvalues, buffer, 1181, 9,
                                   nmax);
    }

    for (size_t m = 0; m < 135; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
