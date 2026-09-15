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


#include "SimdThreeCenterElectronRepulsionGeom010RecGDS.hpp"

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
#include "SimdGeometryI1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdTransferDD.hpp"
#include "SimdTransferDF.hpp"
#include "SimdTransferFD.hpp"
#include "SimdTransferGeom010XDD.hpp"
#include "SimdTransferGeom010XDF.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XFD.hpp"
#include "SimdTransferGeom010XFF.hpp"
#include "SimdTransferGeom010XGD.hpp"
#include "SimdTransferGeom010XPD.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010YDD.hpp"
#include "SimdTransferGeom010YDF.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YFD.hpp"
#include "SimdTransferGeom010YFF.hpp"
#include "SimdTransferGeom010YGD.hpp"
#include "SimdTransferGeom010YPD.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010ZDD.hpp"
#include "SimdTransferGeom010ZDF.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZFD.hpp"
#include "SimdTransferGeom010ZFF.hpp"
#include "SimdTransferGeom010ZGD.hpp"
#include "SimdTransferGeom010ZPD.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_gds_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_gds_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 3261, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 135 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 3261, 577, 556, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 7, ncols,
                                                             fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 15, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 18, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 21, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 24, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 36, 0, 3, 7, 8,
                                                                       15, 18, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 42, 0, 3, 8, 9,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 48, 0, 3, 9, 10,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 54, 0, 3, 10, 11,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 60, 0, 3, 11, 12,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 66, 0, 3, 12, 13,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 72, 0, 3, 15, 18,
                                                                       36, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 82, 0, 3, 18, 21,
                                                                       42, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 92, 0, 3, 21, 24,
                                                                       48, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 102, 0, 3, 24, 27,
                                                                       54, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 112, 0, 3, 27, 30,
                                                                       60, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 122, 0, 3, 36, 42,
                                                                       72, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 137, 0, 3, 42, 48,
                                                                       82, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 152, 0, 3, 48, 54,
                                                                       92, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 167, 0, 3, 54, 60,
                                                                       102, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 182, 0, 3, 72, 82,
                                                                       122, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 203, 0, 3, 82, 92,
                                                                       137, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 224, 0, 3, 92,
                                                                       102, 152, 167, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 245, 0, 3, 122,
                                                                       137, 182, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 273, 0, 3, 137,
                                                                       152, 203, 224, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 301, 0, 3, 182,
                                                                       203, 245, 273, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_d_x(buffer, 337, 15, 72, 1, 1, ncols, beta);

                    simdgeo::geom_d_y(buffer, 343, 15, 72, 1, 1, ncols, beta);

                    simdgeo::geom_d_z(buffer, 349, 15, 72, 1, 1, ncols, beta);

                    simdgeo::geom_f_x(buffer, 355, 36, 122, 1, 1, ncols, beta);

                    simdgeo::geom_f_y(buffer, 365, 36, 122, 1, 1, ncols, beta);

                    simdgeo::geom_f_z(buffer, 375, 36, 122, 1, 1, ncols, beta);

                    simdgeo::geom_g_x(buffer, 385, 72, 182, 1, 1, ncols, beta);

                    simdgeo::geom_g_y(buffer, 400, 72, 182, 1, 1, ncols, beta);

                    simdgeo::geom_g_z(buffer, 415, 72, 182, 1, 1, ncols, beta);

                    simdgeo::geom_h_x(buffer, 430, 122, 245, 1, 1, ncols, beta);

                    simdgeo::geom_h_y(buffer, 451, 122, 245, 1, 1, ncols, beta);

                    simdgeo::geom_h_z(buffer, 472, 122, 245, 1, 1, ncols, beta);

                    simdgeo::geom_i_x(buffer, 493, 182, 301, 1, 1, ncols, beta);

                    simdgeo::geom_i_y(buffer, 521, 182, 301, 1, 1, ncols, beta);

                    simdgeo::geom_i_z(buffer, 549, 182, 301, 1, 1, ncols, beta);

                    simdfunc::contract_primitives(buffer, 577, 337, 6, ncols);

                    simdfunc::contract_primitives(buffer, 589, 343, 6, ncols);

                    simdfunc::contract_primitives(buffer, 601, 349, 6, ncols);

                    simdfunc::contract_primitives(buffer, 613, 36, 6, ncols);

                    simdfunc::contract_primitives(buffer, 625, 355, 10, ncols);

                    simdfunc::contract_primitives(buffer, 645, 365, 10, ncols);

                    simdfunc::contract_primitives(buffer, 665, 375, 10, ncols);

                    simdfunc::contract_primitives(buffer, 685, 72, 10, ncols);

                    simdfunc::contract_primitives(buffer, 705, 385, 15, ncols);

                    simdfunc::contract_primitives(buffer, 735, 400, 15, ncols);

                    simdfunc::contract_primitives(buffer, 765, 415, 15, ncols);

                    simdfunc::contract_primitives(buffer, 795, 122, 15, ncols);

                    simdfunc::contract_primitives(buffer, 825, 430, 21, ncols);

                    simdfunc::contract_primitives(buffer, 867, 451, 21, ncols);

                    simdfunc::contract_primitives(buffer, 909, 472, 21, ncols);

                    simdfunc::contract_primitives(buffer, 951, 182, 21, ncols);

                    simdfunc::contract_primitives(buffer, 993, 493, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1049, 521, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1105, 549, 28, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 583, 577, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 595, 589, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 607, 601, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 619, 613, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 635, 625, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 655, 645, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 675, 665, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 695, 685, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 720, 705, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 750, 735, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 780, 765, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 810, 795, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 846, 825, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 888, 867, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 930, 909, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 972, 951, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1021, 993, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1077, 1049, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1133, 1105, 28, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pd(buffer, coordinates, 1161, 583, 619, 635, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pd(buffer, coordinates, 1179, 595, 619, 655, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pd(buffer, coordinates, 1197, 607, 619, 675, 1, nmax);

        simdtrf::compute_hrr_pd(buffer, coordinates, 1215, 619, 695, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 1233, 635, 695, 720, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 1263, 655, 695, 750, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 1293, 675, 695, 780, 1, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 1323, 695, 810, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 1353, 720, 810, 846, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 1398, 750, 810, 888, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 1443, 780, 810, 930, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 1488, 810, 972, 1, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 1533, 846, 972, 1021, 1, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 1596, 888, 972, 1077, 1, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 1659, 930, 972, 1133, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dd(buffer, coordinates, 1722, 1161, 1215, 1233, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dd(buffer, coordinates, 1758, 1179, 1215, 1263, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dd(buffer, coordinates, 1794, 1197, 1215, 1293, 1, nmax);

        simdtrf::compute_hrr_dd(buffer, coordinates, 1830, 1215, 1323, 1, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 1866, 1233, 1323, 1353, 1, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 1926, 1263, 1323, 1398, 1, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 1986, 1293, 1323, 1443, 1, nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 2046, 1323, 1488, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 2106, 1353, 1488, 1533, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 2196, 1398, 1488, 1596, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 2286, 1443, 1488, 1659, 1, nmax);

        simdtrf::compute_hrr_geom_010x_fd_out_of_second(buffer, coordinates, 2376, 1722, 1830,
                                                        1866, 1, nmax);

        simdtrf::compute_hrr_geom_010y_fd_out_of_second(buffer, coordinates, 2436, 1758, 1830,
                                                        1926, 1, nmax);

        simdtrf::compute_hrr_geom_010z_fd_out_of_second(buffer, coordinates, 2496, 1794, 1830,
                                                        1986, 1, nmax);

        simdtrf::compute_hrr_fd_out_of_second(buffer, coordinates, 2556, 1830, 2046, 1, nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 2616, 1866, 2046, 2106, 1, nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 2716, 1926, 2046, 2196, 1, nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 2816, 1986, 2046, 2286, 1, nmax);

        simdtrf::compute_hrr_geom_010x_gd_out_of_second(buffer, coordinates, 2916, 2376, 2556,
                                                        2616, 1, nmax);

        simdtrf::compute_hrr_geom_010y_gd_out_of_second(buffer, coordinates, 3006, 2436, 2556,
                                                        2716, 1, nmax);

        simdtrf::compute_hrr_geom_010z_gd_out_of_second(buffer, coordinates, 3096, 2496, 2556,
                                                        2816, 1, nmax);

        simdtrf::transform_d_inner(buffer, 3186, 2916, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 3186, 5, nmax);

        simdtrf::transform_d_inner(buffer, 3186, 3006, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 45 * nvalues + n * npairs, nvalues, buffer, 3186, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 3186, 3096, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 90 * nvalues + n * npairs, nvalues, buffer, 3186, 5,
                                   nmax);
    }

    for (size_t m = 0; m < 135; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
