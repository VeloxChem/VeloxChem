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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecDDS.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferGeom100XDD.hpp"
#include "SimdTransferGeom100XDP.hpp"
#include "SimdTransferGeom100XFP.hpp"
#include "SimdTransferGeom100YDD.hpp"
#include "SimdTransferGeom100YDP.hpp"
#include "SimdTransferGeom100YFP.hpp"
#include "SimdTransferGeom100ZDD.hpp"
#include "SimdTransferGeom100ZDP.hpp"
#include "SimdTransferGeom100ZFP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_dds_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_dds_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 1452, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 150 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 1452, 446, 421, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto alpha = a_exps[i];

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fa = -b_exps[j] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pa(buffer, coordinates, 0, nmax, fa);

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

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 20, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 23, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 7, 8,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 8, 9,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 9, 10,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 10, 11,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 14, 15,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 15, 16,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 16, 17,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 17, 18,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 20, 23,
                                                                       50, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 23, 26,
                                                                       56, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 26, 29,
                                                                       62, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 35, 38,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 38, 41,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 41, 44,
                                                                       86, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 50, 56,
                                                                       98, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 56, 62,
                                                                       108, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 74, 80,
                                                                       128, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 80, 86,
                                                                       138, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 98,
                                                                       108, 158, 173, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 239, 0, 3, 128,
                                                                       138, 188, 203, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_d_x(buffer, 260, 20, 98, 1, 1, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 266, 20, 98, 1, 1, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 272, 20, 98, 1, 1, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 278, 35, 128, 1, 1, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 284, 35, 128, 1, 1, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 290, 35, 128, 1, 1, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 296, 50, 158, 1, 1, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 306, 50, 158, 1, 1, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 316, 50, 158, 1, 1, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 326, 74, 188, 1, 1, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 336, 74, 188, 1, 1, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 346, 74, 188, 1, 1, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 356, 98, 218, 1, 1, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 371, 98, 218, 1, 1, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 386, 98, 218, 1, 1, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 401, 128, 239, 1, 1, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 416, 128, 239, 1, 1, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 431, 128, 239, 1, 1, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 446, 260, 6, ncols);

                    simdfunc::contract_primitives(buffer, 458, 266, 6, ncols);

                    simdfunc::contract_primitives(buffer, 470, 272, 6, ncols);

                    simdfunc::contract_primitives(buffer, 482, 50, 6, ncols);

                    simdfunc::contract_primitives(buffer, 494, 278, 6, ncols);

                    simdfunc::contract_primitives(buffer, 506, 284, 6, ncols);

                    simdfunc::contract_primitives(buffer, 518, 290, 6, ncols);

                    simdfunc::contract_primitives(buffer, 530, 74, 6, ncols);

                    simdfunc::contract_primitives(buffer, 542, 296, 10, ncols);

                    simdfunc::contract_primitives(buffer, 562, 306, 10, ncols);

                    simdfunc::contract_primitives(buffer, 582, 316, 10, ncols);

                    simdfunc::contract_primitives(buffer, 602, 98, 10, ncols);

                    simdfunc::contract_primitives(buffer, 622, 326, 10, ncols);

                    simdfunc::contract_primitives(buffer, 642, 336, 10, ncols);

                    simdfunc::contract_primitives(buffer, 662, 346, 10, ncols);

                    simdfunc::contract_primitives(buffer, 682, 128, 10, ncols);

                    simdfunc::contract_primitives(buffer, 702, 356, 15, ncols);

                    simdfunc::contract_primitives(buffer, 732, 371, 15, ncols);

                    simdfunc::contract_primitives(buffer, 762, 386, 15, ncols);

                    simdfunc::contract_primitives(buffer, 792, 401, 15, ncols);

                    simdfunc::contract_primitives(buffer, 822, 416, 15, ncols);

                    simdfunc::contract_primitives(buffer, 852, 431, 15, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 452, 446, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 464, 458, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 476, 470, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 488, 482, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 500, 494, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 512, 506, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 524, 518, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 536, 530, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 552, 542, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 572, 562, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 592, 582, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 612, 602, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 632, 622, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 652, 642, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 672, 662, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 692, 682, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 717, 702, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 747, 732, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 777, 762, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 807, 792, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 837, 822, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 867, 852, 15, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 882, 452, 488, 552,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 900, 464, 488, 572,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 918, 476, 488, 592,
                                                       1, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 936, 488, 612, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 954, 500, 536, 632,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 972, 512, 536, 652,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 990, 524, 536, 672,
                                                       1, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 1008, 536, 692, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 1026, 552, 612, 717,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 1056, 572, 612, 747,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 1086, 592, 612, 777,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 1116, 632, 692, 807,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 1146, 652, 692, 837,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 1176, 672, 692, 867,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 1206, 882, 936, 1026,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 1242, 900, 936, 1056,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 1278, 918, 936, 1086,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 1314, 954, 1008,
                                                       1116, 1, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 1350, 972, 1008,
                                                       1146, 1, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 1386, 990, 1008,
                                                       1176, 1, nmax);

        simdtrf::transform_d_inner(buffer, 1422, 1314, 6, 1, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 1422, 5, nmax);

        simdtrf::transform_d_inner(buffer, 1422, 1350, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 25 * nvalues + n * npairs, nvalues, buffer, 1422, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 1422, 1386, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 50 * nvalues + n * npairs, nvalues, buffer, 1422, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 1422, 1206, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 75 * nvalues + n * npairs, nvalues, buffer, 1422, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 1422, 1242, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 100 * nvalues + n * npairs, nvalues, buffer, 1422, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 1422, 1278, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 125 * nvalues + n * npairs, nvalues, buffer, 1422, 5,
                                   nmax);
    }

    for (size_t m = 0; m < 150; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
