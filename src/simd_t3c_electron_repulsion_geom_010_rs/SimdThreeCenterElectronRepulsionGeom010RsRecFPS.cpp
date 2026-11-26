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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecFPS.hpp"

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
#include "SimdGeometryP1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferGeom010XDD.hpp"
#include "SimdTransferGeom010XDP.hpp"
#include "SimdTransferGeom010XFP.hpp"
#include "SimdTransferGeom010XPD.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPP.hpp"
#include "SimdTransferGeom010YDD.hpp"
#include "SimdTransferGeom010YDP.hpp"
#include "SimdTransferGeom010YFP.hpp"
#include "SimdTransferGeom010YPD.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPP.hpp"
#include "SimdTransferGeom010ZDD.hpp"
#include "SimdTransferGeom010ZDP.hpp"
#include "SimdTransferGeom010ZFP.hpp"
#include "SimdTransferGeom010ZPD.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPP.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransferPP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_fps_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_fps_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 1914, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 1914, 464, 469, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 5,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 13, 3, 5,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 20, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 50, 0, 3, 7, 8,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 56, 0, 3, 8, 9,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 62, 0, 3, 9, 10,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 68, 0, 3, 10, 11,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 74, 0, 3, 14, 15,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 80, 0, 3, 15, 16,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 86, 0, 3, 16, 17,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 92, 0, 3, 17, 18,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 98, 0, 3, 20, 23,
                                                                       50, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 108, 0, 3, 23, 26,
                                                                       56, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 118, 0, 3, 26, 29,
                                                                       62, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 128, 0, 3, 35, 38,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 138, 0, 3, 38, 41,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 148, 0, 3, 41, 44,
                                                                       86, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 158, 0, 3, 50, 56,
                                                                       98, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 173, 0, 3, 56, 62,
                                                                       108, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 188, 0, 3, 74, 80,
                                                                       128, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 80, 86,
                                                                       138, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 218, 0, 3, 98,
                                                                       108, 158, 173, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 239, 0, 3, 128,
                                                                       138, 188, 203, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_p_x(buffer, 260, 7, 50, 1, 1, ncols, beta);

                    simdgeo::geom_p_y(buffer, 263, 7, 50, 1, 1, ncols, beta);

                    simdgeo::geom_p_z(buffer, 266, 7, 50, 1, 1, ncols, beta);

                    simdgeo::geom_p_x(buffer, 269, 14, 74, 1, 1, ncols, beta);

                    simdgeo::geom_p_y(buffer, 272, 14, 74, 1, 1, ncols, beta);

                    simdgeo::geom_p_z(buffer, 275, 14, 74, 1, 1, ncols, beta);

                    simdgeo::geom_d_x(buffer, 278, 20, 98, 1, 1, ncols, beta);

                    simdgeo::geom_d_y(buffer, 284, 20, 98, 1, 1, ncols, beta);

                    simdgeo::geom_d_z(buffer, 290, 20, 98, 1, 1, ncols, beta);

                    simdgeo::geom_d_x(buffer, 296, 35, 128, 1, 1, ncols, beta);

                    simdgeo::geom_d_y(buffer, 302, 35, 128, 1, 1, ncols, beta);

                    simdgeo::geom_d_z(buffer, 308, 35, 128, 1, 1, ncols, beta);

                    simdgeo::geom_f_x(buffer, 314, 50, 158, 1, 1, ncols, beta);

                    simdgeo::geom_f_y(buffer, 324, 50, 158, 1, 1, ncols, beta);

                    simdgeo::geom_f_z(buffer, 334, 50, 158, 1, 1, ncols, beta);

                    simdgeo::geom_f_x(buffer, 344, 74, 188, 1, 1, ncols, beta);

                    simdgeo::geom_f_y(buffer, 354, 74, 188, 1, 1, ncols, beta);

                    simdgeo::geom_f_z(buffer, 364, 74, 188, 1, 1, ncols, beta);

                    simdgeo::geom_g_x(buffer, 374, 98, 218, 1, 1, ncols, beta);

                    simdgeo::geom_g_y(buffer, 389, 98, 218, 1, 1, ncols, beta);

                    simdgeo::geom_g_z(buffer, 404, 98, 218, 1, 1, ncols, beta);

                    simdgeo::geom_g_x(buffer, 419, 128, 239, 1, 1, ncols, beta);

                    simdgeo::geom_g_y(buffer, 434, 128, 239, 1, 1, ncols, beta);

                    simdgeo::geom_g_z(buffer, 449, 128, 239, 1, 1, ncols, beta);

                    simdfunc::contract_primitives(buffer, 464, 260, 3, ncols);

                    simdfunc::contract_primitives(buffer, 470, 263, 3, ncols);

                    simdfunc::contract_primitives(buffer, 476, 266, 3, ncols);

                    simdfunc::contract_primitives(buffer, 482, 20, 3, ncols);

                    simdfunc::contract_primitives(buffer, 488, 269, 3, ncols);

                    simdfunc::contract_primitives(buffer, 494, 272, 3, ncols);

                    simdfunc::contract_primitives(buffer, 500, 275, 3, ncols);

                    simdfunc::contract_primitives(buffer, 506, 35, 3, ncols);

                    simdfunc::contract_primitives(buffer, 512, 278, 6, ncols);

                    simdfunc::contract_primitives(buffer, 524, 284, 6, ncols);

                    simdfunc::contract_primitives(buffer, 536, 290, 6, ncols);

                    simdfunc::contract_primitives(buffer, 548, 50, 6, ncols);

                    simdfunc::contract_primitives(buffer, 560, 296, 6, ncols);

                    simdfunc::contract_primitives(buffer, 572, 302, 6, ncols);

                    simdfunc::contract_primitives(buffer, 584, 308, 6, ncols);

                    simdfunc::contract_primitives(buffer, 596, 74, 6, ncols);

                    simdfunc::contract_primitives(buffer, 608, 314, 10, ncols);

                    simdfunc::contract_primitives(buffer, 628, 324, 10, ncols);

                    simdfunc::contract_primitives(buffer, 648, 334, 10, ncols);

                    simdfunc::contract_primitives(buffer, 668, 98, 10, ncols);

                    simdfunc::contract_primitives(buffer, 688, 344, 10, ncols);

                    simdfunc::contract_primitives(buffer, 708, 354, 10, ncols);

                    simdfunc::contract_primitives(buffer, 728, 364, 10, ncols);

                    simdfunc::contract_primitives(buffer, 748, 128, 10, ncols);

                    simdfunc::contract_primitives(buffer, 768, 374, 15, ncols);

                    simdfunc::contract_primitives(buffer, 798, 389, 15, ncols);

                    simdfunc::contract_primitives(buffer, 828, 404, 15, ncols);

                    simdfunc::contract_primitives(buffer, 858, 419, 15, ncols);

                    simdfunc::contract_primitives(buffer, 888, 434, 15, ncols);

                    simdfunc::contract_primitives(buffer, 918, 449, 15, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 467, 464, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 473, 470, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 479, 476, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 485, 482, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 491, 488, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 497, 494, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 503, 500, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 509, 506, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 518, 512, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 530, 524, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 542, 536, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 554, 548, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 566, 560, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 578, 572, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 590, 584, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 602, 596, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 618, 608, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 638, 628, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 658, 648, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 678, 668, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 698, 688, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 718, 708, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 738, 728, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 758, 748, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 783, 768, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 813, 798, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 843, 828, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 873, 858, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 903, 888, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 933, 918, 15, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 948, 467, 485, 518, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 957, 473, 485, 530, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 966, 479, 485, 542, 1, nmax);

        simdtrf::compute_hrr_pp(buffer, coordinates, 975, 485, 554, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 984, 491, 509, 566, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 993, 497, 509, 578, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 1002, 503, 509, 590, 1, nmax);

        simdtrf::compute_hrr_pp(buffer, coordinates, 1011, 509, 602, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 1020, 518, 554, 618, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 1038, 530, 554, 638, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 1056, 542, 554, 658, 1, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 1074, 554, 678, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 1092, 566, 602, 698, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 1110, 578, 602, 718, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 1128, 590, 602, 738, 1, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 1146, 602, 758, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 1164, 618, 678, 783, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 1194, 638, 678, 813, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 1224, 658, 678, 843, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 1254, 698, 758, 873, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 1284, 718, 758, 903, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 1314, 738, 758, 933, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dp_out_of_second(buffer, coordinates, 1344, 948, 975,
                                                        1020, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dp_out_of_second(buffer, coordinates, 1362, 957, 975,
                                                        1038, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dp_out_of_second(buffer, coordinates, 1380, 966, 975,
                                                        1056, 1, nmax);

        simdtrf::compute_hrr_dp_out_of_second(buffer, coordinates, 1398, 975, 1074, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dp_out_of_second(buffer, coordinates, 1416, 984, 1011,
                                                        1092, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dp_out_of_second(buffer, coordinates, 1434, 993, 1011,
                                                        1110, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dp_out_of_second(buffer, coordinates, 1452, 1002, 1011,
                                                        1128, 1, nmax);

        simdtrf::compute_hrr_dp_out_of_second(buffer, coordinates, 1470, 1011, 1146, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 1488, 1020, 1074, 1164, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 1524, 1038, 1074, 1194, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 1560, 1056, 1074, 1224, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 1596, 1092, 1146, 1254, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 1632, 1110, 1146, 1284, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 1668, 1128, 1146, 1314, 1, nmax);

        simdtrf::compute_hrr_geom_010x_fp_out_of_second(buffer, coordinates, 1704, 1344, 1398,
                                                        1488, 1, nmax);

        simdtrf::compute_hrr_geom_010y_fp_out_of_second(buffer, coordinates, 1734, 1362, 1398,
                                                        1524, 1, nmax);

        simdtrf::compute_hrr_geom_010z_fp_out_of_second(buffer, coordinates, 1764, 1380, 1398,
                                                        1560, 1, nmax);

        simdtrf::compute_hrr_geom_010x_fp_out_of_second(buffer, coordinates, 1794, 1416, 1470,
                                                        1596, 1, nmax);

        simdtrf::compute_hrr_geom_010y_fp_out_of_second(buffer, coordinates, 1824, 1434, 1470,
                                                        1632, 1, nmax);

        simdtrf::compute_hrr_geom_010z_fp_out_of_second(buffer, coordinates, 1854, 1452, 1470,
                                                        1668, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1884, 1794, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 1884, 3, nmax);

        simdtrf::transform_p_inner(buffer, 1884, 1824, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 21 * nvalues + n * npairs, nvalues, buffer, 1884, 3,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 1884, 1854, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 42 * nvalues + n * npairs, nvalues, buffer, 1884, 3,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 1884, 1704, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 1884, 3,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 1884, 1734, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 84 * nvalues + n * npairs, nvalues, buffer, 1884, 3,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 1884, 1764, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 105 * nvalues + n * npairs, nvalues, buffer, 1884, 3,
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
