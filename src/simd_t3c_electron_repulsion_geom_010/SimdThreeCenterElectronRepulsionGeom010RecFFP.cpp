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


#include "SimdThreeCenterElectronRepulsionGeom010RecFFP.hpp"

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
#include "SimdGeometryI1.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSKP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSKS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdTransferDF.hpp"
#include "SimdTransferGeom010XDF.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XFF.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010YDF.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YFF.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010ZDF.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZFF.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_ffp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_ffp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 7066, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 441 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 7066, 1351, 1524, dimensions);

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
                                                        5, 6, 7, 8}, ncols, fj, i * nprim_b + j,
                                                        fq);

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

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 337, 3, 15, 36,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 355, 3, 36, 72,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 385, 3, 72, 122,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 430, 3, 122, 182,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 493, 3, 182, 245,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 577, 3, 245, 301,
                                                                       ncols, p, q);

                    simdgeo::geom_f_x(buffer, 685, 337, 385, 1, 3, ncols, beta);

                    simdgeo::geom_f_y(buffer, 715, 337, 385, 1, 3, ncols, beta);

                    simdgeo::geom_f_z(buffer, 745, 337, 385, 1, 3, ncols, beta);

                    simdgeo::geom_g_x(buffer, 775, 355, 430, 1, 3, ncols, beta);

                    simdgeo::geom_g_y(buffer, 820, 355, 430, 1, 3, ncols, beta);

                    simdgeo::geom_g_z(buffer, 865, 355, 430, 1, 3, ncols, beta);

                    simdgeo::geom_h_x(buffer, 910, 385, 493, 1, 3, ncols, beta);

                    simdgeo::geom_h_y(buffer, 973, 385, 493, 1, 3, ncols, beta);

                    simdgeo::geom_h_z(buffer, 1036, 385, 493, 1, 3, ncols, beta);

                    simdgeo::geom_i_x(buffer, 1099, 430, 577, 1, 3, ncols, beta);

                    simdgeo::geom_i_y(buffer, 1183, 430, 577, 1, 3, ncols, beta);

                    simdgeo::geom_i_z(buffer, 1267, 430, 577, 1, 3, ncols, beta);

                    simdfunc::contract_primitives(buffer, 1351, 685, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1411, 715, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1471, 745, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1531, 355, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1591, 775, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1681, 820, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1771, 865, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1861, 385, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1951, 910, 63, ncols);

                    simdfunc::contract_primitives(buffer, 2077, 973, 63, ncols);

                    simdfunc::contract_primitives(buffer, 2203, 1036, 63, ncols);

                    simdfunc::contract_primitives(buffer, 2329, 430, 63, ncols);

                    simdfunc::contract_primitives(buffer, 2455, 1099, 84, ncols);

                    simdfunc::contract_primitives(buffer, 2623, 1183, 84, ncols);

                    simdfunc::contract_primitives(buffer, 2791, 1267, 84, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 1381, 1351, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1441, 1411, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1501, 1471, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1561, 1531, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1636, 1591, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1726, 1681, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1816, 1771, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1906, 1861, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2014, 1951, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2140, 2077, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2266, 2203, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2392, 2329, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2539, 2455, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2707, 2623, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2875, 2791, 28, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 2959, 1381, 1561, 1636, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 3049, 1441, 1561, 1726, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 3139, 1501, 1561, 1816, 3, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 3229, 1561, 1906, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 3319, 1636, 1906, 2014, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 3454, 1726, 1906, 2140, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 3589, 1816, 1906, 2266, 3, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 3724, 1906, 2392, 3, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 3859, 2014, 2392, 2539, 3, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 4048, 2140, 2392, 2707, 3, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 4237, 2266, 2392, 2875, 3, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 4426, 2959, 3229, 3319, 3, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 4606, 3049, 3229, 3454, 3, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 4786, 3139, 3229, 3589, 3, nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 4966, 3229, 3724, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 5146, 3319, 3724, 3859, 3, nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 5416, 3454, 3724, 4048, 3, nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 5686, 3589, 3724, 4237, 3, nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 5956, 4426, 4966, 5146, 3, nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 6256, 4606, 4966, 5416, 3, nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 6556, 4786, 4966, 5686, 3, nmax);

        simdtrf::transform_f_inner(buffer, 6856, 5956, 10, 3, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 6856, 21, nmax);

        simdtrf::transform_f_inner(buffer, 6856, 6256, 10, 3, nmax);

        simdtrf::transform_f_outer(values + 147 * nvalues + n * npairs, nvalues, buffer, 6856,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 6856, 6556, 10, 3, nmax);

        simdtrf::transform_f_outer(values + 294 * nvalues + n * npairs, nvalues, buffer, 6856,
                                   21, nmax);
    }

    for (size_t m = 0; m < 441; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
