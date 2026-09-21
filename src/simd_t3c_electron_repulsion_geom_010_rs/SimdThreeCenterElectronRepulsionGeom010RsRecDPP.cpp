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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecDPP.hpp"

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
compute_rs_geom_010_dpp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_dpp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2410, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 270 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 2410, 700, 762, dimensions);

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
                                                            4, 5}, ncols, fj, i * nprim_b + j,
                                                            fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 12, 3, {1, 2, 3, 4,
                                                        5}, ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 18, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 42, 0, 3, 7, 8,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 48, 0, 3, 8, 9,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 54, 0, 3, 9, 10,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 60, 0, 3, 13, 14,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 66, 0, 3, 14, 15,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 72, 0, 3, 15, 16,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 78, 0, 3, 18, 21,
                                                                       42, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 88, 0, 3, 21, 24,
                                                                       48, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 98, 0, 3, 30, 33,
                                                                       60, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 108, 0, 3, 33, 36,
                                                                       66, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 118, 0, 3, 42, 48,
                                                                       78, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 133, 0, 3, 60, 66,
                                                                       98, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 148, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 151, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 154, 3, 7, 18,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 163, 3, 13, 30,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 172, 3, 18, 42,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 190, 3, 30, 60,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 208, 3, 42, 78,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 238, 3, 60, 98,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 268, 3, 78, 118,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 313, 3, 98, 133,
                                                                       ncols, p, q);

                    simdgeo::geom_p_x(buffer, 358, 148, 172, 1, 3, ncols, beta);

                    simdgeo::geom_p_y(buffer, 367, 148, 172, 1, 3, ncols, beta);

                    simdgeo::geom_p_z(buffer, 376, 148, 172, 1, 3, ncols, beta);

                    simdgeo::geom_p_x(buffer, 385, 151, 190, 1, 3, ncols, beta);

                    simdgeo::geom_p_y(buffer, 394, 151, 190, 1, 3, ncols, beta);

                    simdgeo::geom_p_z(buffer, 403, 151, 190, 1, 3, ncols, beta);

                    simdgeo::geom_d_x(buffer, 412, 154, 208, 1, 3, ncols, beta);

                    simdgeo::geom_d_y(buffer, 430, 154, 208, 1, 3, ncols, beta);

                    simdgeo::geom_d_z(buffer, 448, 154, 208, 1, 3, ncols, beta);

                    simdgeo::geom_d_x(buffer, 466, 163, 238, 1, 3, ncols, beta);

                    simdgeo::geom_d_y(buffer, 484, 163, 238, 1, 3, ncols, beta);

                    simdgeo::geom_d_z(buffer, 502, 163, 238, 1, 3, ncols, beta);

                    simdgeo::geom_f_x(buffer, 520, 172, 268, 1, 3, ncols, beta);

                    simdgeo::geom_f_y(buffer, 550, 172, 268, 1, 3, ncols, beta);

                    simdgeo::geom_f_z(buffer, 580, 172, 268, 1, 3, ncols, beta);

                    simdgeo::geom_f_x(buffer, 610, 190, 313, 1, 3, ncols, beta);

                    simdgeo::geom_f_y(buffer, 640, 190, 313, 1, 3, ncols, beta);

                    simdgeo::geom_f_z(buffer, 670, 190, 313, 1, 3, ncols, beta);

                    simdfunc::contract_primitives(buffer, 700, 358, 9, ncols);

                    simdfunc::contract_primitives(buffer, 718, 367, 9, ncols);

                    simdfunc::contract_primitives(buffer, 736, 376, 9, ncols);

                    simdfunc::contract_primitives(buffer, 754, 154, 9, ncols);

                    simdfunc::contract_primitives(buffer, 772, 385, 9, ncols);

                    simdfunc::contract_primitives(buffer, 790, 394, 9, ncols);

                    simdfunc::contract_primitives(buffer, 808, 403, 9, ncols);

                    simdfunc::contract_primitives(buffer, 826, 163, 9, ncols);

                    simdfunc::contract_primitives(buffer, 844, 412, 18, ncols);

                    simdfunc::contract_primitives(buffer, 880, 430, 18, ncols);

                    simdfunc::contract_primitives(buffer, 916, 448, 18, ncols);

                    simdfunc::contract_primitives(buffer, 952, 172, 18, ncols);

                    simdfunc::contract_primitives(buffer, 988, 466, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1024, 484, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1060, 502, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1096, 190, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1132, 520, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1192, 550, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1252, 580, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1312, 610, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1372, 640, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1432, 670, 30, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 709, 700, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 727, 718, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 745, 736, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 763, 754, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 781, 772, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 799, 790, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 817, 808, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 835, 826, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 862, 844, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 898, 880, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 934, 916, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 970, 952, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1006, 988, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1042, 1024, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1078, 1060, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1114, 1096, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1162, 1132, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1222, 1192, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1282, 1252, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1342, 1312, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1402, 1372, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1462, 1432, 10, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 1492, 709, 763, 862, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 1519, 727, 763, 898, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 1546, 745, 763, 934, 3, nmax);

        simdtrf::compute_hrr_pp(buffer, coordinates, 1573, 763, 970, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pp(buffer, coordinates, 1600, 781, 835, 1006, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pp(buffer, coordinates, 1627, 799, 835, 1042, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pp(buffer, coordinates, 1654, 817, 835, 1078, 3, nmax);

        simdtrf::compute_hrr_pp(buffer, coordinates, 1681, 835, 1114, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 1708, 862, 970, 1162, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 1762, 898, 970, 1222, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 1816, 934, 970, 1282, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 1870, 1006, 1114, 1342, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 1924, 1042, 1114, 1402, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 1978, 1078, 1114, 1462, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dp_out_of_second(buffer, coordinates, 2032, 1492, 1573,
                                                        1708, 3, nmax);

        simdtrf::compute_hrr_geom_010y_dp_out_of_second(buffer, coordinates, 2086, 1519, 1573,
                                                        1762, 3, nmax);

        simdtrf::compute_hrr_geom_010z_dp_out_of_second(buffer, coordinates, 2140, 1546, 1573,
                                                        1816, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dp_out_of_second(buffer, coordinates, 2194, 1600, 1681,
                                                        1870, 3, nmax);

        simdtrf::compute_hrr_geom_010y_dp_out_of_second(buffer, coordinates, 2248, 1627, 1681,
                                                        1924, 3, nmax);

        simdtrf::compute_hrr_geom_010z_dp_out_of_second(buffer, coordinates, 2302, 1654, 1681,
                                                        1978, 3, nmax);

        simdtrf::transform_p_inner(buffer, 2356, 2194, 6, 3, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 2356, 9, nmax);

        simdtrf::transform_p_inner(buffer, 2356, 2248, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 45 * nvalues + n * npairs, nvalues, buffer, 2356, 9,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 2356, 2302, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 90 * nvalues + n * npairs, nvalues, buffer, 2356, 9,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 2356, 2032, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 135 * nvalues + n * npairs, nvalues, buffer, 2356, 9,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 2356, 2086, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 180 * nvalues + n * npairs, nvalues, buffer, 2356, 9,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 2356, 2140, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 225 * nvalues + n * npairs, nvalues, buffer, 2356, 9,
                                   nmax);
    }

    for (size_t m = 0; m < 270; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
