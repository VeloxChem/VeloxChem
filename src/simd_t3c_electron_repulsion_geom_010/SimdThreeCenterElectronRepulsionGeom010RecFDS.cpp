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


#include "SimdThreeCenterElectronRepulsionGeom010RecFDS.hpp"

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
#include "SimdGeometryG1.hpp"
#include "SimdGeometryH1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdTransferDD.hpp"
#include "SimdTransferGeom010XDD.hpp"
#include "SimdTransferGeom010XDF.hpp"
#include "SimdTransferGeom010XFD.hpp"
#include "SimdTransferGeom010XPD.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010YDD.hpp"
#include "SimdTransferGeom010YDF.hpp"
#include "SimdTransferGeom010YFD.hpp"
#include "SimdTransferGeom010YPD.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010ZDD.hpp"
#include "SimdTransferGeom010ZDF.hpp"
#include "SimdTransferGeom010ZFD.hpp"
#include "SimdTransferGeom010ZPD.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_fds_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_fds_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 1628, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 105 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 1628, 373, 353, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 6, ncols,
                                                             fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 14, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 17, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 20, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 32, 0, 3, 7, 8,
                                                                       14, 17, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 38, 0, 3, 8, 9,
                                                                       17, 20, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 44, 0, 3, 9, 10,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 50, 0, 3, 10, 11,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 56, 0, 3, 11, 12,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 62, 0, 3, 14, 17,
                                                                       32, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 72, 0, 3, 17, 20,
                                                                       38, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 82, 0, 3, 20, 23,
                                                                       44, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 92, 0, 3, 23, 26,
                                                                       50, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 102, 0, 3, 32, 38,
                                                                       62, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 117, 0, 3, 38, 44,
                                                                       72, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 132, 0, 3, 44, 50,
                                                                       82, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 147, 0, 3, 62, 72,
                                                                       102, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 168, 0, 3, 72, 82,
                                                                       117, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 189, 0, 3, 102,
                                                                       117, 147, 168, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_d_x(buffer, 217, 14, 62, 1, 1, ncols, beta);

                    simdgeo::geom_d_y(buffer, 223, 14, 62, 1, 1, ncols, beta);

                    simdgeo::geom_d_z(buffer, 229, 14, 62, 1, 1, ncols, beta);

                    simdgeo::geom_f_x(buffer, 235, 32, 102, 1, 1, ncols, beta);

                    simdgeo::geom_f_y(buffer, 245, 32, 102, 1, 1, ncols, beta);

                    simdgeo::geom_f_z(buffer, 255, 32, 102, 1, 1, ncols, beta);

                    simdgeo::geom_g_x(buffer, 265, 62, 147, 1, 1, ncols, beta);

                    simdgeo::geom_g_y(buffer, 280, 62, 147, 1, 1, ncols, beta);

                    simdgeo::geom_g_z(buffer, 295, 62, 147, 1, 1, ncols, beta);

                    simdgeo::geom_h_x(buffer, 310, 102, 189, 1, 1, ncols, beta);

                    simdgeo::geom_h_y(buffer, 331, 102, 189, 1, 1, ncols, beta);

                    simdgeo::geom_h_z(buffer, 352, 102, 189, 1, 1, ncols, beta);

                    simdfunc::contract_primitives(buffer, 373, 217, 6, ncols);

                    simdfunc::contract_primitives(buffer, 385, 223, 6, ncols);

                    simdfunc::contract_primitives(buffer, 397, 229, 6, ncols);

                    simdfunc::contract_primitives(buffer, 409, 32, 6, ncols);

                    simdfunc::contract_primitives(buffer, 421, 235, 10, ncols);

                    simdfunc::contract_primitives(buffer, 441, 245, 10, ncols);

                    simdfunc::contract_primitives(buffer, 461, 255, 10, ncols);

                    simdfunc::contract_primitives(buffer, 481, 62, 10, ncols);

                    simdfunc::contract_primitives(buffer, 501, 265, 15, ncols);

                    simdfunc::contract_primitives(buffer, 531, 280, 15, ncols);

                    simdfunc::contract_primitives(buffer, 561, 295, 15, ncols);

                    simdfunc::contract_primitives(buffer, 591, 102, 15, ncols);

                    simdfunc::contract_primitives(buffer, 621, 310, 21, ncols);

                    simdfunc::contract_primitives(buffer, 663, 331, 21, ncols);

                    simdfunc::contract_primitives(buffer, 705, 352, 21, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 379, 373, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 391, 385, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 403, 397, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 415, 409, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 431, 421, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 451, 441, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 471, 461, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 491, 481, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 516, 501, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 546, 531, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 576, 561, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 606, 591, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 642, 621, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 684, 663, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 726, 705, 21, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 747, 379, 415, 431, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 765, 391, 415, 451, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 783, 403, 415, 471, 1, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 801, 415, 491, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 819, 431, 491, 516, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 849, 451, 491, 546, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 879, 471, 491, 576, 1, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 909, 491, 606, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 939, 516, 606, 642, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 984, 546, 606, 684, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 1029, 576, 606, 726, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 1074, 747, 801, 819, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 1110, 765, 801, 849, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 1146, 783, 801, 879, 1, nmax);

        simdtrf::compute_hrr_dd(buffer, coordinates, 1182, 801, 909, 1, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 1218, 819, 909, 939, 1, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 1278, 849, 909, 984, 1, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 1338, 879, 909, 1029, 1, nmax);

        simdtrf::compute_hrr_geom_010x_fd_out_of_second(buffer, coordinates, 1398, 1074, 1182,
                                                        1218, 1, nmax);

        simdtrf::compute_hrr_geom_010y_fd_out_of_second(buffer, coordinates, 1458, 1110, 1182,
                                                        1278, 1, nmax);

        simdtrf::compute_hrr_geom_010z_fd_out_of_second(buffer, coordinates, 1518, 1146, 1182,
                                                        1338, 1, nmax);

        simdtrf::transform_d_inner(buffer, 1578, 1398, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 1578, 5, nmax);

        simdtrf::transform_d_inner(buffer, 1578, 1458, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 35 * nvalues + n * npairs, nvalues, buffer, 1578, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 1578, 1518, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 70 * nvalues + n * npairs, nvalues, buffer, 1578, 5,
                                   nmax);
    }

    for (size_t m = 0; m < 105; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
