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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecDFS.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdGeometryF1.hpp"
#include "SimdGeometryG1.hpp"
#include "SimdGeometryH1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdTransferGeom010XDF.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010YDF.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010ZDF.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_dfs_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_dfs_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2268, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 210 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 2268, 704, 631, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 6,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 14, 3, 6,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 58, 0, 3, 7, 8,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 64, 0, 3, 8, 9,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 70, 0, 3, 9, 10,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 76, 0, 3, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 82, 0, 3, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 88, 0, 3, 15, 16,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 94, 0, 3, 16, 17,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 100, 0, 3, 17, 18,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 106, 0, 3, 18, 19,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 112, 0, 3, 19, 20,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 118, 0, 3, 22, 25,
                                                                       58, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 128, 0, 3, 25, 28,
                                                                       64, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 138, 0, 3, 28, 31,
                                                                       70, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 148, 0, 3, 31, 34,
                                                                       76, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 158, 0, 3, 40, 43,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 168, 0, 3, 43, 46,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 178, 0, 3, 46, 49,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 188, 0, 3, 49, 52,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 198, 0, 3, 58, 64,
                                                                       118, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 213, 0, 3, 64, 70,
                                                                       128, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 228, 0, 3, 70, 76,
                                                                       138, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 243, 0, 3, 88, 94,
                                                                       158, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 258, 0, 3, 94,
                                                                       100, 168, 178, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 100,
                                                                       106, 178, 188, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 288, 0, 3, 118,
                                                                       128, 198, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 309, 0, 3, 128,
                                                                       138, 213, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 330, 0, 3, 158,
                                                                       168, 243, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 351, 0, 3, 168,
                                                                       178, 258, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 372, 0, 3, 198,
                                                                       213, 288, 309, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 400, 0, 3, 243,
                                                                       258, 330, 351, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_f_x(buffer, 428, 58, 198, 1, 1, ncols, beta);

                    simdgeo::geom_f_y(buffer, 438, 58, 198, 1, 1, ncols, beta);

                    simdgeo::geom_f_z(buffer, 448, 58, 198, 1, 1, ncols, beta);

                    simdgeo::geom_f_x(buffer, 458, 88, 243, 1, 1, ncols, beta);

                    simdgeo::geom_f_y(buffer, 468, 88, 243, 1, 1, ncols, beta);

                    simdgeo::geom_f_z(buffer, 478, 88, 243, 1, 1, ncols, beta);

                    simdgeo::geom_g_x(buffer, 488, 118, 288, 1, 1, ncols, beta);

                    simdgeo::geom_g_y(buffer, 503, 118, 288, 1, 1, ncols, beta);

                    simdgeo::geom_g_z(buffer, 518, 118, 288, 1, 1, ncols, beta);

                    simdgeo::geom_g_x(buffer, 533, 158, 330, 1, 1, ncols, beta);

                    simdgeo::geom_g_y(buffer, 548, 158, 330, 1, 1, ncols, beta);

                    simdgeo::geom_g_z(buffer, 563, 158, 330, 1, 1, ncols, beta);

                    simdgeo::geom_h_x(buffer, 578, 198, 372, 1, 1, ncols, beta);

                    simdgeo::geom_h_y(buffer, 599, 198, 372, 1, 1, ncols, beta);

                    simdgeo::geom_h_z(buffer, 620, 198, 372, 1, 1, ncols, beta);

                    simdgeo::geom_h_x(buffer, 641, 243, 400, 1, 1, ncols, beta);

                    simdgeo::geom_h_y(buffer, 662, 243, 400, 1, 1, ncols, beta);

                    simdgeo::geom_h_z(buffer, 683, 243, 400, 1, 1, ncols, beta);

                    simdfunc::contract_primitives(buffer, 704, 428, 10, ncols);

                    simdfunc::contract_primitives(buffer, 724, 438, 10, ncols);

                    simdfunc::contract_primitives(buffer, 744, 448, 10, ncols);

                    simdfunc::contract_primitives(buffer, 764, 118, 10, ncols);

                    simdfunc::contract_primitives(buffer, 784, 458, 10, ncols);

                    simdfunc::contract_primitives(buffer, 804, 468, 10, ncols);

                    simdfunc::contract_primitives(buffer, 824, 478, 10, ncols);

                    simdfunc::contract_primitives(buffer, 844, 158, 10, ncols);

                    simdfunc::contract_primitives(buffer, 864, 488, 15, ncols);

                    simdfunc::contract_primitives(buffer, 894, 503, 15, ncols);

                    simdfunc::contract_primitives(buffer, 924, 518, 15, ncols);

                    simdfunc::contract_primitives(buffer, 954, 198, 15, ncols);

                    simdfunc::contract_primitives(buffer, 984, 533, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1014, 548, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1044, 563, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1074, 243, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1104, 578, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1146, 599, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1188, 620, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1230, 641, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1272, 662, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1314, 683, 21, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 714, 704, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 734, 724, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 754, 744, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 774, 764, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 794, 784, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 814, 804, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 834, 824, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 854, 844, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 879, 864, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 909, 894, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 939, 924, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 969, 954, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 999, 984, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1029, 1014, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1059, 1044, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1089, 1074, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1125, 1104, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1167, 1146, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1209, 1188, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1251, 1230, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1293, 1272, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1335, 1314, 21, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 1356, 714, 774, 879, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 1386, 734, 774, 909, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 1416, 754, 774, 939, 1, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 1446, 774, 969, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 1476, 794, 854, 999, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 1506, 814, 854, 1029, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 1536, 834, 854, 1059, 1, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 1566, 854, 1089, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 1596, 879, 969, 1125, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 1641, 909, 969, 1167, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 1686, 939, 969, 1209, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 1731, 999, 1089, 1251, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 1776, 1029, 1089, 1293, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 1821, 1059, 1089, 1335, 1, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 1866, 1356, 1446, 1596, 1, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 1926, 1386, 1446, 1641, 1, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 1986, 1416, 1446, 1686, 1, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 2046, 1476, 1566, 1731, 1, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 2106, 1506, 1566, 1776, 1, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 2166, 1536, 1566, 1821, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2226, 2046, 6, 1, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 2226, 7, nmax);

        simdtrf::transform_f_inner(buffer, 2226, 2106, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 35 * nvalues + n * npairs, nvalues, buffer, 2226, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 2226, 2166, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 70 * nvalues + n * npairs, nvalues, buffer, 2226, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 2226, 1866, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 105 * nvalues + n * npairs, nvalues, buffer, 2226, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 2226, 1926, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 140 * nvalues + n * npairs, nvalues, buffer, 2226, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 2226, 1986, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 175 * nvalues + n * npairs, nvalues, buffer, 2226, 7,
                                   nmax);
    }

    for (size_t m = 0; m < 210; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
