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


#include "SimdThreeCenterElectronRepulsionGeom010RecFGS.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdGeometryG1.hpp"
#include "SimdGeometryH1.hpp"
#include "SimdGeometryI1.hpp"
#include "SimdGeometryK1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdTransferDG.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XDH.hpp"
#include "SimdTransferGeom010XFG.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010XPI.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YDH.hpp"
#include "SimdTransferGeom010YFG.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010YPI.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZDH.hpp"
#include "SimdTransferGeom010ZFG.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferGeom010ZPI.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_fgs_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_fgs_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 3492, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 189 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 3492, 802, 692, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 8, ncols,
                                                             fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 16, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 19, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 40, 0, 3, 7, 8,
                                                                       16, 19, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 46, 0, 3, 8, 9,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 52, 0, 3, 9, 10,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 58, 0, 3, 10, 11,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 64, 0, 3, 11, 12,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 70, 0, 3, 12, 13,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 76, 0, 3, 13, 14,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 82, 0, 3, 16, 19,
                                                                       40, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 92, 0, 3, 19, 22,
                                                                       46, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 102, 0, 3, 22, 25,
                                                                       52, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 112, 0, 3, 25, 28,
                                                                       58, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 122, 0, 3, 28, 31,
                                                                       64, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 132, 0, 3, 31, 34,
                                                                       70, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 142, 0, 3, 40, 46,
                                                                       82, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 157, 0, 3, 46, 52,
                                                                       92, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 172, 0, 3, 52, 58,
                                                                       102, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 187, 0, 3, 58, 64,
                                                                       112, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 202, 0, 3, 64, 70,
                                                                       122, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 217, 0, 3, 82, 92,
                                                                       142, 157, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 238, 0, 3, 92,
                                                                       102, 157, 172, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 259, 0, 3, 102,
                                                                       112, 172, 187, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 280, 0, 3, 112,
                                                                       122, 187, 202, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 301, 0, 3, 142,
                                                                       157, 217, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 329, 0, 3, 157,
                                                                       172, 238, 259, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 357, 0, 3, 172,
                                                                       187, 259, 280, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 385, 0, 3, 217,
                                                                       238, 301, 329, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 421, 0, 3, 238,
                                                                       259, 329, 357, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 457, 0, 3, 301,
                                                                       329, 385, 421, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_g_x(buffer, 502, 82, 217, 1, 1, ncols, beta);

                    simdgeo::geom_g_y(buffer, 517, 82, 217, 1, 1, ncols, beta);

                    simdgeo::geom_g_z(buffer, 532, 82, 217, 1, 1, ncols, beta);

                    simdgeo::geom_h_x(buffer, 547, 142, 301, 1, 1, ncols, beta);

                    simdgeo::geom_h_y(buffer, 568, 142, 301, 1, 1, ncols, beta);

                    simdgeo::geom_h_z(buffer, 589, 142, 301, 1, 1, ncols, beta);

                    simdgeo::geom_i_x(buffer, 610, 217, 385, 1, 1, ncols, beta);

                    simdgeo::geom_i_y(buffer, 638, 217, 385, 1, 1, ncols, beta);

                    simdgeo::geom_i_z(buffer, 666, 217, 385, 1, 1, ncols, beta);

                    simdgeo::geom_k_x(buffer, 694, 301, 457, 1, 1, ncols, beta);

                    simdgeo::geom_k_y(buffer, 730, 301, 457, 1, 1, ncols, beta);

                    simdgeo::geom_k_z(buffer, 766, 301, 457, 1, 1, ncols, beta);

                    simdfunc::contract_primitives(buffer, 802, 502, 15, ncols);

                    simdfunc::contract_primitives(buffer, 832, 517, 15, ncols);

                    simdfunc::contract_primitives(buffer, 862, 532, 15, ncols);

                    simdfunc::contract_primitives(buffer, 892, 142, 15, ncols);

                    simdfunc::contract_primitives(buffer, 922, 547, 21, ncols);

                    simdfunc::contract_primitives(buffer, 964, 568, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1006, 589, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1048, 217, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1090, 610, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1146, 638, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1202, 666, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1258, 301, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1314, 694, 36, ncols);

                    simdfunc::contract_primitives(buffer, 1386, 730, 36, ncols);

                    simdfunc::contract_primitives(buffer, 1458, 766, 36, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 817, 802, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 847, 832, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 877, 862, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 907, 892, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 943, 922, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 985, 964, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1027, 1006, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1069, 1048, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1118, 1090, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1174, 1146, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1230, 1202, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1286, 1258, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1350, 1314, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1422, 1386, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1494, 1458, 36, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 1530, 817, 907, 943, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 1575, 847, 907, 985, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 1620, 877, 907, 1027, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 1665, 907, 1069, 1, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 1710, 943, 1069, 1118, 1, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 1773, 985, 1069, 1174, 1, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 1836, 1027, 1069, 1230, 1, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 1899, 1069, 1286, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 1962, 1118, 1286, 1350, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 2046, 1174, 1286, 1422, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 2130, 1230, 1286, 1494, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 2214, 1530, 1665, 1710, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 2304, 1575, 1665, 1773, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 2394, 1620, 1665, 1836, 1, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 2484, 1665, 1899, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 2574, 1710, 1899, 1962, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 2700, 1773, 1899, 2046, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 2826, 1836, 1899, 2130, 1, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 2952, 2214, 2484, 2574, 1, nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 3102, 2304, 2484, 2700, 1, nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 3252, 2394, 2484, 2826, 1, nmax);

        simdtrf::transform_g_inner(buffer, 3402, 2952, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 3402, 9, nmax);

        simdtrf::transform_g_inner(buffer, 3402, 3102, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 3402, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 3402, 3252, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 126 * nvalues + n * npairs, nvalues, buffer, 3402, 9,
                                   nmax);
    }

    for (size_t m = 0; m < 189; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
