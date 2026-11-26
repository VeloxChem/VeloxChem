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


#include "SimdThreeCenterElectronRepulsionGeom010RecDFP.hpp"

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
compute_geom_010_dfp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_dfp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 3280, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 315 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 3280, 871, 915, dimensions);

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
                                                        5, 6, 7}, ncols, fj, i * nprim_b + j,
                                                        fq);

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

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 217, 3, 14, 32,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 235, 3, 32, 62,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 265, 3, 62, 102,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 310, 3, 102, 147,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 373, 3, 147, 189,
                                                                       ncols, p, q);

                    simdgeo::geom_f_x(buffer, 457, 217, 265, 1, 3, ncols, beta);

                    simdgeo::geom_f_y(buffer, 487, 217, 265, 1, 3, ncols, beta);

                    simdgeo::geom_f_z(buffer, 517, 217, 265, 1, 3, ncols, beta);

                    simdgeo::geom_g_x(buffer, 547, 235, 310, 1, 3, ncols, beta);

                    simdgeo::geom_g_y(buffer, 592, 235, 310, 1, 3, ncols, beta);

                    simdgeo::geom_g_z(buffer, 637, 235, 310, 1, 3, ncols, beta);

                    simdgeo::geom_h_x(buffer, 682, 265, 373, 1, 3, ncols, beta);

                    simdgeo::geom_h_y(buffer, 745, 265, 373, 1, 3, ncols, beta);

                    simdgeo::geom_h_z(buffer, 808, 265, 373, 1, 3, ncols, beta);

                    simdfunc::contract_primitives(buffer, 871, 457, 30, ncols);

                    simdfunc::contract_primitives(buffer, 931, 487, 30, ncols);

                    simdfunc::contract_primitives(buffer, 991, 517, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1051, 235, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1111, 547, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1201, 592, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1291, 637, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1381, 265, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1471, 682, 63, ncols);

                    simdfunc::contract_primitives(buffer, 1597, 745, 63, ncols);

                    simdfunc::contract_primitives(buffer, 1723, 808, 63, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 901, 871, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 961, 931, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1021, 991, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1081, 1051, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1156, 1111, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1246, 1201, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1336, 1291, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1426, 1381, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1534, 1471, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1660, 1597, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1786, 1723, 21, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 1849, 901, 1081, 1156, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 1939, 961, 1081, 1246, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 2029, 1021, 1081, 1336, 3, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 2119, 1081, 1426, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 2209, 1156, 1426, 1534, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 2344, 1246, 1426, 1660, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 2479, 1336, 1426, 1786, 3, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 2614, 1849, 2119, 2209, 3, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 2794, 1939, 2119, 2344, 3, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 2974, 2029, 2119, 2479, 3, nmax);

        simdtrf::transform_f_inner(buffer, 3154, 2614, 6, 3, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 3154, 21, nmax);

        simdtrf::transform_f_inner(buffer, 3154, 2794, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 105 * nvalues + n * npairs, nvalues, buffer, 3154,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 3154, 2974, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 210 * nvalues + n * npairs, nvalues, buffer, 3154,
                                   21, nmax);
    }

    for (size_t m = 0; m < 315; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
