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


#include "SimdThreeCenterElectronRepulsionGeom010RecGFP.hpp"

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
#include "SimdGeometryK1.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSLP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSLS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdTransferDF.hpp"
#include "SimdTransferDG.hpp"
#include "SimdTransferFF.hpp"
#include "SimdTransferGeom010XDF.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XDH.hpp"
#include "SimdTransferGeom010XFF.hpp"
#include "SimdTransferGeom010XFG.hpp"
#include "SimdTransferGeom010XGF.hpp"
#include "SimdTransferGeom010XPF.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010XPI.hpp"
#include "SimdTransferGeom010YDF.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YDH.hpp"
#include "SimdTransferGeom010YFF.hpp"
#include "SimdTransferGeom010YFG.hpp"
#include "SimdTransferGeom010YGF.hpp"
#include "SimdTransferGeom010YPF.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010YPI.hpp"
#include "SimdTransferGeom010ZDF.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZDH.hpp"
#include "SimdTransferGeom010ZFF.hpp"
#include "SimdTransferGeom010ZFG.hpp"
#include "SimdTransferGeom010ZGF.hpp"
#include "SimdTransferGeom010ZPF.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferGeom010ZPI.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_gfp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_gfp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 13960, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 567 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 13960, 1975, 2316, dimensions);

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
                                                        5, 6, 7, 8, 9}, ncols, fj,
                                                        i * nprim_b + j, fq);

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

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 502, 3, 16, 40,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 520, 3, 40, 82,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 550, 3, 82, 142,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 595, 3, 142, 217,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 658, 3, 217, 301,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 742, 3, 301, 385,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 850, 3, 385, 457,
                                                                       ncols, p, q);

                    simdgeo::geom_f_x(buffer, 985, 502, 550, 1, 3, ncols, beta);

                    simdgeo::geom_f_y(buffer, 1015, 502, 550, 1, 3, ncols, beta);

                    simdgeo::geom_f_z(buffer, 1045, 502, 550, 1, 3, ncols, beta);

                    simdgeo::geom_g_x(buffer, 1075, 520, 595, 1, 3, ncols, beta);

                    simdgeo::geom_g_y(buffer, 1120, 520, 595, 1, 3, ncols, beta);

                    simdgeo::geom_g_z(buffer, 1165, 520, 595, 1, 3, ncols, beta);

                    simdgeo::geom_h_x(buffer, 1210, 550, 658, 1, 3, ncols, beta);

                    simdgeo::geom_h_y(buffer, 1273, 550, 658, 1, 3, ncols, beta);

                    simdgeo::geom_h_z(buffer, 1336, 550, 658, 1, 3, ncols, beta);

                    simdgeo::geom_i_x(buffer, 1399, 595, 742, 1, 3, ncols, beta);

                    simdgeo::geom_i_y(buffer, 1483, 595, 742, 1, 3, ncols, beta);

                    simdgeo::geom_i_z(buffer, 1567, 595, 742, 1, 3, ncols, beta);

                    simdgeo::geom_k_x(buffer, 1651, 658, 850, 1, 3, ncols, beta);

                    simdgeo::geom_k_y(buffer, 1759, 658, 850, 1, 3, ncols, beta);

                    simdgeo::geom_k_z(buffer, 1867, 658, 850, 1, 3, ncols, beta);

                    simdfunc::contract_primitives(buffer, 1975, 985, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2035, 1015, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2095, 1045, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2155, 520, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2215, 1075, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2305, 1120, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2395, 1165, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2485, 550, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2575, 1210, 63, ncols);

                    simdfunc::contract_primitives(buffer, 2701, 1273, 63, ncols);

                    simdfunc::contract_primitives(buffer, 2827, 1336, 63, ncols);

                    simdfunc::contract_primitives(buffer, 2953, 595, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3079, 1399, 84, ncols);

                    simdfunc::contract_primitives(buffer, 3247, 1483, 84, ncols);

                    simdfunc::contract_primitives(buffer, 3415, 1567, 84, ncols);

                    simdfunc::contract_primitives(buffer, 3583, 658, 84, ncols);

                    simdfunc::contract_primitives(buffer, 3751, 1651, 108, ncols);

                    simdfunc::contract_primitives(buffer, 3967, 1759, 108, ncols);

                    simdfunc::contract_primitives(buffer, 4183, 1867, 108, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 2005, 1975, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2065, 2035, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2125, 2095, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2185, 2155, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2260, 2215, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2350, 2305, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2440, 2395, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2530, 2485, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2638, 2575, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2764, 2701, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2890, 2827, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3016, 2953, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3163, 3079, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3331, 3247, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3499, 3415, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3667, 3583, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3859, 3751, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4075, 3967, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4291, 4183, 36, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pf(buffer, coordinates, 4399, 2005, 2185, 2260, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pf(buffer, coordinates, 4489, 2065, 2185, 2350, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pf(buffer, coordinates, 4579, 2125, 2185, 2440, 3, nmax);

        simdtrf::compute_hrr_pf(buffer, coordinates, 4669, 2185, 2530, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 4759, 2260, 2530, 2638, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 4894, 2350, 2530, 2764, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 5029, 2440, 2530, 2890, 3, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 5164, 2530, 3016, 3, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 5299, 2638, 3016, 3163, 3, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 5488, 2764, 3016, 3331, 3, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 5677, 2890, 3016, 3499, 3, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 5866, 3016, 3667, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 6055, 3163, 3667, 3859, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 6307, 3331, 3667, 4075, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 6559, 3499, 3667, 4291, 3, nmax);

        simdtrf::compute_hrr_geom_010x_df(buffer, coordinates, 6811, 4399, 4669, 4759, 3, nmax);

        simdtrf::compute_hrr_geom_010y_df(buffer, coordinates, 6991, 4489, 4669, 4894, 3, nmax);

        simdtrf::compute_hrr_geom_010z_df(buffer, coordinates, 7171, 4579, 4669, 5029, 3, nmax);

        simdtrf::compute_hrr_df(buffer, coordinates, 7351, 4669, 5164, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 7531, 4759, 5164, 5299, 3, nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 7801, 4894, 5164, 5488, 3, nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 8071, 5029, 5164, 5677, 3, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 8341, 5164, 5866, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 8611, 5299, 5866, 6055, 3, nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 8989, 5488, 5866, 6307, 3, nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 9367, 5677, 5866, 6559, 3, nmax);

        simdtrf::compute_hrr_geom_010x_ff(buffer, coordinates, 9745, 6811, 7351, 7531, 3, nmax);

        simdtrf::compute_hrr_geom_010y_ff(buffer, coordinates, 10045, 6991, 7351, 7801, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_ff(buffer, coordinates, 10345, 7171, 7351, 8071, 3,
                                          nmax);

        simdtrf::compute_hrr_ff(buffer, coordinates, 10645, 7351, 8341, 3, nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 10945, 7531, 8341, 8611, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 11395, 7801, 8341, 8989, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 11845, 8071, 8341, 9367, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gf_out_of_second(buffer, coordinates, 12295, 9745, 10645,
                                                        10945, 3, nmax);

        simdtrf::compute_hrr_geom_010y_gf_out_of_second(buffer, coordinates, 12745, 10045, 10645,
                                                        11395, 3, nmax);

        simdtrf::compute_hrr_geom_010z_gf_out_of_second(buffer, coordinates, 13195, 10345, 10645,
                                                        11845, 3, nmax);

        simdtrf::transform_f_inner(buffer, 13645, 12295, 15, 3, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 13645, 21, nmax);

        simdtrf::transform_f_inner(buffer, 13645, 12745, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 189 * nvalues + n * npairs, nvalues, buffer, 13645,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 13645, 13195, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 378 * nvalues + n * npairs, nvalues, buffer, 13645,
                                   21, nmax);
    }

    for (size_t m = 0; m < 567; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
