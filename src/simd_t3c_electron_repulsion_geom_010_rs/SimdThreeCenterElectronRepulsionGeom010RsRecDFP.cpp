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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecDFP.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIP.hpp"
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
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_dfp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_dfp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 6428, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 630 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 6428, 1736, 1893, dimensions);

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
                                                            4, 5, 6, 7}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 14, 3, {1, 2, 3, 4,
                                                        5, 6, 7}, ncols, fj, i * nprim_b + j,
                                                        fq);

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

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 428, 3, 22, 58,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 446, 3, 40, 88,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 464, 3, 58, 118,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 494, 3, 88, 158,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 524, 3, 118, 198,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 569, 3, 158, 243,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 614, 3, 198, 288,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 677, 3, 243, 330,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 740, 3, 288, 372,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 824, 3, 330, 400,
                                                                       ncols, p, q);

                    simdgeo::geom_f_x(buffer, 908, 428, 524, 1, 3, ncols, beta);

                    simdgeo::geom_f_y(buffer, 938, 428, 524, 1, 3, ncols, beta);

                    simdgeo::geom_f_z(buffer, 968, 428, 524, 1, 3, ncols, beta);

                    simdgeo::geom_f_x(buffer, 998, 446, 569, 1, 3, ncols, beta);

                    simdgeo::geom_f_y(buffer, 1028, 446, 569, 1, 3, ncols, beta);

                    simdgeo::geom_f_z(buffer, 1058, 446, 569, 1, 3, ncols, beta);

                    simdgeo::geom_g_x(buffer, 1088, 464, 614, 1, 3, ncols, beta);

                    simdgeo::geom_g_y(buffer, 1133, 464, 614, 1, 3, ncols, beta);

                    simdgeo::geom_g_z(buffer, 1178, 464, 614, 1, 3, ncols, beta);

                    simdgeo::geom_g_x(buffer, 1223, 494, 677, 1, 3, ncols, beta);

                    simdgeo::geom_g_y(buffer, 1268, 494, 677, 1, 3, ncols, beta);

                    simdgeo::geom_g_z(buffer, 1313, 494, 677, 1, 3, ncols, beta);

                    simdgeo::geom_h_x(buffer, 1358, 524, 740, 1, 3, ncols, beta);

                    simdgeo::geom_h_y(buffer, 1421, 524, 740, 1, 3, ncols, beta);

                    simdgeo::geom_h_z(buffer, 1484, 524, 740, 1, 3, ncols, beta);

                    simdgeo::geom_h_x(buffer, 1547, 569, 824, 1, 3, ncols, beta);

                    simdgeo::geom_h_y(buffer, 1610, 569, 824, 1, 3, ncols, beta);

                    simdgeo::geom_h_z(buffer, 1673, 569, 824, 1, 3, ncols, beta);

                    simdfunc::contract_primitives(buffer, 1736, 908, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1796, 938, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1856, 968, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1916, 464, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1976, 998, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2036, 1028, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2096, 1058, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2156, 494, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2216, 1088, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2306, 1133, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2396, 1178, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2486, 524, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2576, 1223, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2666, 1268, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2756, 1313, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2846, 569, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2936, 1358, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3062, 1421, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3188, 1484, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3314, 1547, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3440, 1610, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3566, 1673, 63, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 1766, 1736, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1826, 1796, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1886, 1856, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1946, 1916, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2006, 1976, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2066, 2036, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2126, 2096, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2186, 2156, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2261, 2216, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2351, 2306, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2441, 2396, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2531, 2486, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2621, 2576, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2711, 2666, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2801, 2756, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2891, 2846, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2999, 2936, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3125, 3062, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3251, 3188, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3377, 3314, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3503, 3440, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3629, 3566, 21, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 3692, 1766, 1946, 2261, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 3782, 1826, 1946, 2351, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 3872, 1886, 1946, 2441, 3, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 3962, 1946, 2531, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 4052, 2006, 2186, 2621, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 4142, 2066, 2186, 2711, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 4232, 2126, 2186, 2801, 3, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 4322, 2186, 2891, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 4412, 2261, 2531, 2999, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 4547, 2351, 2531, 3125, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 4682, 2441, 2531, 3251, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 4817, 2621, 2891, 3377, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 4952, 2711, 2891, 3503, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 5087, 2801, 2891, 3629, 3, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 5222, 3692, 3962, 4412, 3, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 5402, 3782, 3962, 4547, 3, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 5582, 3872, 3962, 4682, 3, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 5762, 4052, 4322, 4817, 3, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 5942, 4142, 4322, 4952, 3, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 6122, 4232, 4322, 5087, 3, nmax);

        simdtrf::transform_f_inner(buffer, 6302, 5762, 6, 3, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 6302, 21, nmax);

        simdtrf::transform_f_inner(buffer, 6302, 5942, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 105 * nvalues + n * npairs, nvalues, buffer, 6302,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 6302, 6122, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 210 * nvalues + n * npairs, nvalues, buffer, 6302,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 6302, 5222, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 6302,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 6302, 5402, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 420 * nvalues + n * npairs, nvalues, buffer, 6302,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 6302, 5582, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 525 * nvalues + n * npairs, nvalues, buffer, 6302,
                                   21, nmax);
    }

    for (size_t m = 0; m < 630; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
