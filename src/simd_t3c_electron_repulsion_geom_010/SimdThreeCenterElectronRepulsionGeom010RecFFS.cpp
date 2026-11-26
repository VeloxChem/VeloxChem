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


#include "SimdThreeCenterElectronRepulsionGeom010RecFFS.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSHS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSIS.hpp"
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
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_ffs_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_ffs_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2464, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 147 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 2464, 559, 508, dimensions);

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

                    simdgeo::geom_f_x(buffer, 337, 36, 122, 1, 1, ncols, beta);

                    simdgeo::geom_f_y(buffer, 347, 36, 122, 1, 1, ncols, beta);

                    simdgeo::geom_f_z(buffer, 357, 36, 122, 1, 1, ncols, beta);

                    simdgeo::geom_g_x(buffer, 367, 72, 182, 1, 1, ncols, beta);

                    simdgeo::geom_g_y(buffer, 382, 72, 182, 1, 1, ncols, beta);

                    simdgeo::geom_g_z(buffer, 397, 72, 182, 1, 1, ncols, beta);

                    simdgeo::geom_h_x(buffer, 412, 122, 245, 1, 1, ncols, beta);

                    simdgeo::geom_h_y(buffer, 433, 122, 245, 1, 1, ncols, beta);

                    simdgeo::geom_h_z(buffer, 454, 122, 245, 1, 1, ncols, beta);

                    simdgeo::geom_i_x(buffer, 475, 182, 301, 1, 1, ncols, beta);

                    simdgeo::geom_i_y(buffer, 503, 182, 301, 1, 1, ncols, beta);

                    simdgeo::geom_i_z(buffer, 531, 182, 301, 1, 1, ncols, beta);

                    simdfunc::contract_primitives(buffer, 559, 337, 10, ncols);

                    simdfunc::contract_primitives(buffer, 579, 347, 10, ncols);

                    simdfunc::contract_primitives(buffer, 599, 357, 10, ncols);

                    simdfunc::contract_primitives(buffer, 619, 72, 10, ncols);

                    simdfunc::contract_primitives(buffer, 639, 367, 15, ncols);

                    simdfunc::contract_primitives(buffer, 669, 382, 15, ncols);

                    simdfunc::contract_primitives(buffer, 699, 397, 15, ncols);

                    simdfunc::contract_primitives(buffer, 729, 122, 15, ncols);

                    simdfunc::contract_primitives(buffer, 759, 412, 21, ncols);

                    simdfunc::contract_primitives(buffer, 801, 433, 21, ncols);

                    simdfunc::contract_primitives(buffer, 843, 454, 21, ncols);

                    simdfunc::contract_primitives(buffer, 885, 182, 21, ncols);

                    simdfunc::contract_primitives(buffer, 927, 475, 28, ncols);

                    simdfunc::contract_primitives(buffer, 983, 503, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1039, 531, 28, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 569, 559, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 589, 579, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 609, 599, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 629, 619, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 654, 639, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 684, 669, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 714, 699, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 744, 729, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 780, 759, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 822, 801, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 864, 843, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 906, 885, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 955, 927, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1011, 983, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1067, 1039, 28, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 1095, 569, 629, 654, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 1125, 589, 629, 684, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 1155, 609, 629, 714, 1, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 1185, 629, 744, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 1215, 654, 744, 780, 1, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 1260, 684, 744, 822, 1, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 1305, 714, 744, 864, 1, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 1350, 744, 906, 1, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 1395, 780, 906, 955, 1, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 1458, 822, 906, 1011, 1, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 1521, 864, 906, 1067, 1, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 1584, 1095, 1185, 1215, 1, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 1644, 1125, 1185, 1260, 1, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 1704, 1155, 1185, 1305, 1, nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 1764, 1185, 1350, 1, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 1824, 1215, 1350, 1395, 1, nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 1914, 1260, 1350, 1458, 1, nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 2004, 1305, 1350, 1521, 1, nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 2094, 1584, 1764, 1824, 1, nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 2194, 1644, 1764, 1914, 1, nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 2294, 1704, 1764, 2004, 1, nmax);

        simdtrf::transform_f_inner(buffer, 2394, 2094, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 2394, 7, nmax);

        simdtrf::transform_f_inner(buffer, 2394, 2194, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 49 * nvalues + n * npairs, nvalues, buffer, 2394, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 2394, 2294, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 98 * nvalues + n * npairs, nvalues, buffer, 2394, 7,
                                   nmax);
    }

    for (size_t m = 0; m < 147; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
