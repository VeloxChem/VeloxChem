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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecFDS.hpp"

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
compute_rs_geom_010_fds_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_fds_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 3200, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 3200, 740, 727, dimensions);

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

                    simdgeo::geom_d_x(buffer, 428, 22, 118, 1, 1, ncols, beta);

                    simdgeo::geom_d_y(buffer, 434, 22, 118, 1, 1, ncols, beta);

                    simdgeo::geom_d_z(buffer, 440, 22, 118, 1, 1, ncols, beta);

                    simdgeo::geom_d_x(buffer, 446, 40, 158, 1, 1, ncols, beta);

                    simdgeo::geom_d_y(buffer, 452, 40, 158, 1, 1, ncols, beta);

                    simdgeo::geom_d_z(buffer, 458, 40, 158, 1, 1, ncols, beta);

                    simdgeo::geom_f_x(buffer, 464, 58, 198, 1, 1, ncols, beta);

                    simdgeo::geom_f_y(buffer, 474, 58, 198, 1, 1, ncols, beta);

                    simdgeo::geom_f_z(buffer, 484, 58, 198, 1, 1, ncols, beta);

                    simdgeo::geom_f_x(buffer, 494, 88, 243, 1, 1, ncols, beta);

                    simdgeo::geom_f_y(buffer, 504, 88, 243, 1, 1, ncols, beta);

                    simdgeo::geom_f_z(buffer, 514, 88, 243, 1, 1, ncols, beta);

                    simdgeo::geom_g_x(buffer, 524, 118, 288, 1, 1, ncols, beta);

                    simdgeo::geom_g_y(buffer, 539, 118, 288, 1, 1, ncols, beta);

                    simdgeo::geom_g_z(buffer, 554, 118, 288, 1, 1, ncols, beta);

                    simdgeo::geom_g_x(buffer, 569, 158, 330, 1, 1, ncols, beta);

                    simdgeo::geom_g_y(buffer, 584, 158, 330, 1, 1, ncols, beta);

                    simdgeo::geom_g_z(buffer, 599, 158, 330, 1, 1, ncols, beta);

                    simdgeo::geom_h_x(buffer, 614, 198, 372, 1, 1, ncols, beta);

                    simdgeo::geom_h_y(buffer, 635, 198, 372, 1, 1, ncols, beta);

                    simdgeo::geom_h_z(buffer, 656, 198, 372, 1, 1, ncols, beta);

                    simdgeo::geom_h_x(buffer, 677, 243, 400, 1, 1, ncols, beta);

                    simdgeo::geom_h_y(buffer, 698, 243, 400, 1, 1, ncols, beta);

                    simdgeo::geom_h_z(buffer, 719, 243, 400, 1, 1, ncols, beta);

                    simdfunc::contract_primitives(buffer, 740, 428, 6, ncols);

                    simdfunc::contract_primitives(buffer, 752, 434, 6, ncols);

                    simdfunc::contract_primitives(buffer, 764, 440, 6, ncols);

                    simdfunc::contract_primitives(buffer, 776, 58, 6, ncols);

                    simdfunc::contract_primitives(buffer, 788, 446, 6, ncols);

                    simdfunc::contract_primitives(buffer, 800, 452, 6, ncols);

                    simdfunc::contract_primitives(buffer, 812, 458, 6, ncols);

                    simdfunc::contract_primitives(buffer, 824, 88, 6, ncols);

                    simdfunc::contract_primitives(buffer, 836, 464, 10, ncols);

                    simdfunc::contract_primitives(buffer, 856, 474, 10, ncols);

                    simdfunc::contract_primitives(buffer, 876, 484, 10, ncols);

                    simdfunc::contract_primitives(buffer, 896, 118, 10, ncols);

                    simdfunc::contract_primitives(buffer, 916, 494, 10, ncols);

                    simdfunc::contract_primitives(buffer, 936, 504, 10, ncols);

                    simdfunc::contract_primitives(buffer, 956, 514, 10, ncols);

                    simdfunc::contract_primitives(buffer, 976, 158, 10, ncols);

                    simdfunc::contract_primitives(buffer, 996, 524, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1026, 539, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1056, 554, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1086, 198, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1116, 569, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1146, 584, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1176, 599, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1206, 243, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1236, 614, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1278, 635, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1320, 656, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1362, 677, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1404, 698, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1446, 719, 21, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 746, 740, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 758, 752, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 770, 764, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 782, 776, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 794, 788, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 806, 800, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 818, 812, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 830, 824, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 846, 836, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 866, 856, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 886, 876, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 906, 896, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 926, 916, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 946, 936, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 966, 956, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 986, 976, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1011, 996, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1041, 1026, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1071, 1056, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1101, 1086, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1131, 1116, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1161, 1146, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1191, 1176, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1221, 1206, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1257, 1236, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1299, 1278, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1341, 1320, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1383, 1362, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1425, 1404, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1467, 1446, 21, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 1488, 746, 782, 846, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 1506, 758, 782, 866, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 1524, 770, 782, 886, 1, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 1542, 782, 906, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 1560, 794, 830, 926, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 1578, 806, 830, 946, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 1596, 818, 830, 966, 1, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 1614, 830, 986, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 1632, 846, 906, 1011, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 1662, 866, 906, 1041, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 1692, 886, 906, 1071, 1, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 1722, 906, 1101, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 1752, 926, 986, 1131, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 1782, 946, 986, 1161, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 1812, 966, 986, 1191, 1, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 1842, 986, 1221, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 1872, 1011, 1101, 1257, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 1917, 1041, 1101, 1299, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 1962, 1071, 1101, 1341, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 2007, 1131, 1221, 1383, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 2052, 1161, 1221, 1425, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 2097, 1191, 1221, 1467, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 2142, 1488, 1542, 1632, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 2178, 1506, 1542, 1662, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 2214, 1524, 1542, 1692, 1, nmax);

        simdtrf::compute_hrr_dd(buffer, coordinates, 2250, 1542, 1722, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 2286, 1560, 1614, 1752, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 2322, 1578, 1614, 1782, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 2358, 1596, 1614, 1812, 1, nmax);

        simdtrf::compute_hrr_dd(buffer, coordinates, 2394, 1614, 1842, 1, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 2430, 1632, 1722, 1872, 1, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 2490, 1662, 1722, 1917, 1, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 2550, 1692, 1722, 1962, 1, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 2610, 1752, 1842, 2007, 1, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 2670, 1782, 1842, 2052, 1, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 2730, 1812, 1842, 2097, 1, nmax);

        simdtrf::compute_hrr_geom_010x_fd_out_of_second(buffer, coordinates, 2790, 2142, 2250,
                                                        2430, 1, nmax);

        simdtrf::compute_hrr_geom_010y_fd_out_of_second(buffer, coordinates, 2850, 2178, 2250,
                                                        2490, 1, nmax);

        simdtrf::compute_hrr_geom_010z_fd_out_of_second(buffer, coordinates, 2910, 2214, 2250,
                                                        2550, 1, nmax);

        simdtrf::compute_hrr_geom_010x_fd_out_of_second(buffer, coordinates, 2970, 2286, 2394,
                                                        2610, 1, nmax);

        simdtrf::compute_hrr_geom_010y_fd_out_of_second(buffer, coordinates, 3030, 2322, 2394,
                                                        2670, 1, nmax);

        simdtrf::compute_hrr_geom_010z_fd_out_of_second(buffer, coordinates, 3090, 2358, 2394,
                                                        2730, 1, nmax);

        simdtrf::transform_d_inner(buffer, 3150, 2970, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 3150, 5, nmax);

        simdtrf::transform_d_inner(buffer, 3150, 3030, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 35 * nvalues + n * npairs, nvalues, buffer, 3150, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 3150, 3090, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 70 * nvalues + n * npairs, nvalues, buffer, 3150, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 3150, 2790, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 105 * nvalues + n * npairs, nvalues, buffer, 3150, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 3150, 2850, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 140 * nvalues + n * npairs, nvalues, buffer, 3150, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 3150, 2910, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 175 * nvalues + n * npairs, nvalues, buffer, 3150, 5,
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
