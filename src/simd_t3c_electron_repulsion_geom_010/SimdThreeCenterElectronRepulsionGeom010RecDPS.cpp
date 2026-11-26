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


#include "SimdThreeCenterElectronRepulsionGeom010RecDPS.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
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
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_dps_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_dps_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 428, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 45 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 428, 134, 122, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 4, ncols,
                                                             fj, i * nprim_b + j, fq);

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

                    simdgeo::geom_p_x(buffer, 77, 7, 24, 1, 1, ncols, beta);

                    simdgeo::geom_p_y(buffer, 80, 7, 24, 1, 1, ncols, beta);

                    simdgeo::geom_p_z(buffer, 83, 7, 24, 1, 1, ncols, beta);

                    simdgeo::geom_d_x(buffer, 86, 12, 42, 1, 1, ncols, beta);

                    simdgeo::geom_d_y(buffer, 92, 12, 42, 1, 1, ncols, beta);

                    simdgeo::geom_d_z(buffer, 98, 12, 42, 1, 1, ncols, beta);

                    simdgeo::geom_f_x(buffer, 104, 24, 62, 1, 1, ncols, beta);

                    simdgeo::geom_f_y(buffer, 114, 24, 62, 1, 1, ncols, beta);

                    simdgeo::geom_f_z(buffer, 124, 24, 62, 1, 1, ncols, beta);

                    simdfunc::contract_primitives(buffer, 134, 77, 3, ncols);

                    simdfunc::contract_primitives(buffer, 140, 80, 3, ncols);

                    simdfunc::contract_primitives(buffer, 146, 83, 3, ncols);

                    simdfunc::contract_primitives(buffer, 152, 12, 3, ncols);

                    simdfunc::contract_primitives(buffer, 158, 86, 6, ncols);

                    simdfunc::contract_primitives(buffer, 170, 92, 6, ncols);

                    simdfunc::contract_primitives(buffer, 182, 98, 6, ncols);

                    simdfunc::contract_primitives(buffer, 194, 24, 6, ncols);

                    simdfunc::contract_primitives(buffer, 206, 104, 10, ncols);

                    simdfunc::contract_primitives(buffer, 226, 114, 10, ncols);

                    simdfunc::contract_primitives(buffer, 246, 124, 10, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 137, 134, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 143, 140, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 149, 146, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 155, 152, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 164, 158, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 176, 170, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 188, 182, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 200, 194, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 216, 206, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 236, 226, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 256, 246, 10, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 266, 137, 155, 164, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 275, 143, 155, 176, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 284, 149, 155, 188, 1, nmax);

        simdtrf::compute_hrr_pp(buffer, coordinates, 293, 155, 200, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 302, 164, 200, 216, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 320, 176, 200, 236, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 338, 188, 200, 256, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dp_out_of_second(buffer, coordinates, 356, 266, 293, 302,
                                                        1, nmax);

        simdtrf::compute_hrr_geom_010y_dp_out_of_second(buffer, coordinates, 374, 275, 293, 320,
                                                        1, nmax);

        simdtrf::compute_hrr_geom_010z_dp_out_of_second(buffer, coordinates, 392, 284, 293, 338,
                                                        1, nmax);

        simdtrf::transform_p_inner(buffer, 410, 356, 6, 1, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 410, 3, nmax);

        simdtrf::transform_p_inner(buffer, 410, 374, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 15 * nvalues + n * npairs, nvalues, buffer, 410, 3,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 410, 392, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 30 * nvalues + n * npairs, nvalues, buffer, 410, 3,
                                   nmax);
    }

    for (size_t m = 0; m < 45; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
