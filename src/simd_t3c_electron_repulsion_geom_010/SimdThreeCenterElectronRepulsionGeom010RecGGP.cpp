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


#include "SimdThreeCenterElectronRepulsionGeom010RecGGP.hpp"

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
#include "SimdGeometryL1.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSMP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSMS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdTransferDG.hpp"
#include "SimdTransferDH.hpp"
#include "SimdTransferFG.hpp"
#include "SimdTransferGeom010XDG.hpp"
#include "SimdTransferGeom010XDH.hpp"
#include "SimdTransferGeom010XDI.hpp"
#include "SimdTransferGeom010XFG.hpp"
#include "SimdTransferGeom010XFH.hpp"
#include "SimdTransferGeom010XGG.hpp"
#include "SimdTransferGeom010XPG.hpp"
#include "SimdTransferGeom010XPH.hpp"
#include "SimdTransferGeom010XPI.hpp"
#include "SimdTransferGeom010XPK.hpp"
#include "SimdTransferGeom010YDG.hpp"
#include "SimdTransferGeom010YDH.hpp"
#include "SimdTransferGeom010YDI.hpp"
#include "SimdTransferGeom010YFG.hpp"
#include "SimdTransferGeom010YFH.hpp"
#include "SimdTransferGeom010YGG.hpp"
#include "SimdTransferGeom010YPG.hpp"
#include "SimdTransferGeom010YPH.hpp"
#include "SimdTransferGeom010YPI.hpp"
#include "SimdTransferGeom010YPK.hpp"
#include "SimdTransferGeom010ZDG.hpp"
#include "SimdTransferGeom010ZDH.hpp"
#include "SimdTransferGeom010ZDI.hpp"
#include "SimdTransferGeom010ZFG.hpp"
#include "SimdTransferGeom010ZFH.hpp"
#include "SimdTransferGeom010ZGG.hpp"
#include "SimdTransferGeom010ZPG.hpp"
#include "SimdTransferGeom010ZPH.hpp"
#include "SimdTransferGeom010ZPI.hpp"
#include "SimdTransferGeom010ZPK.hpp"
#include "SimdTransferPG.hpp"
#include "SimdTransferPH.hpp"
#include "SimdTransferPI.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_010_ggp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_ggp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 19367, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 729 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 19367, 2657, 3075, dimensions);

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
                                                        5, 6, 7, 8, 9, 10}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 17, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 20, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 44, 0, 3, 7, 8,
                                                                       17, 20, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 50, 0, 3, 8, 9,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 56, 0, 3, 9, 10,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 62, 0, 3, 10, 11,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 68, 0, 3, 11, 12,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 74, 0, 3, 12, 13,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 80, 0, 3, 13, 14,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 86, 0, 3, 14, 15,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 92, 0, 3, 17, 20,
                                                                       44, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 102, 0, 3, 20, 23,
                                                                       50, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 112, 0, 3, 23, 26,
                                                                       56, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 122, 0, 3, 26, 29,
                                                                       62, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 132, 0, 3, 29, 32,
                                                                       68, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 142, 0, 3, 32, 35,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 152, 0, 3, 35, 38,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 162, 0, 3, 44, 50,
                                                                       92, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 177, 0, 3, 50, 56,
                                                                       102, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 192, 0, 3, 56, 62,
                                                                       112, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 207, 0, 3, 62, 68,
                                                                       122, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 222, 0, 3, 68, 74,
                                                                       132, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 237, 0, 3, 74, 80,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 252, 0, 3, 92,
                                                                       102, 162, 177, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 273, 0, 3, 102,
                                                                       112, 177, 192, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 294, 0, 3, 112,
                                                                       122, 192, 207, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 315, 0, 3, 122,
                                                                       132, 207, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_shs_three_center_electron_repulsion_0(buffer, 336, 0, 3, 132,
                                                                       142, 222, 237, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 357, 0, 3, 162,
                                                                       177, 252, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 385, 0, 3, 177,
                                                                       192, 273, 294, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 413, 0, 3, 192,
                                                                       207, 294, 315, ncols,
                                                                       gamma, p, q);

                    compute_prim_sis_three_center_electron_repulsion_0(buffer, 441, 0, 3, 207,
                                                                       222, 315, 336, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 469, 0, 3, 252,
                                                                       273, 357, 385, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 505, 0, 3, 273,
                                                                       294, 385, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_sks_three_center_electron_repulsion_0(buffer, 541, 0, 3, 294,
                                                                       315, 413, 441, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 577, 0, 3, 357,
                                                                       385, 469, 505, ncols,
                                                                       gamma, p, q);

                    compute_prim_sls_three_center_electron_repulsion_0(buffer, 622, 0, 3, 385,
                                                                       413, 505, 541, ncols,
                                                                       gamma, p, q);

                    compute_prim_sms_three_center_electron_repulsion_0(buffer, 667, 0, 3, 469,
                                                                       505, 577, 622, ncols,
                                                                       gamma, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 722, 3, 44, 92,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 752, 3, 92, 162,
                                                                       ncols, p, q);

                    compute_prim_shp_three_center_electron_repulsion_0(buffer, 797, 3, 162, 252,
                                                                       ncols, p, q);

                    compute_prim_sip_three_center_electron_repulsion_0(buffer, 860, 3, 252, 357,
                                                                       ncols, p, q);

                    compute_prim_skp_three_center_electron_repulsion_0(buffer, 944, 3, 357, 469,
                                                                       ncols, p, q);

                    compute_prim_slp_three_center_electron_repulsion_0(buffer, 1052, 3, 469, 577,
                                                                       ncols, p, q);

                    compute_prim_smp_three_center_electron_repulsion_0(buffer, 1187, 3, 577, 667,
                                                                       ncols, p, q);

                    simdgeo::geom_g_x(buffer, 1352, 722, 797, 1, 3, ncols, beta);

                    simdgeo::geom_g_y(buffer, 1397, 722, 797, 1, 3, ncols, beta);

                    simdgeo::geom_g_z(buffer, 1442, 722, 797, 1, 3, ncols, beta);

                    simdgeo::geom_h_x(buffer, 1487, 752, 860, 1, 3, ncols, beta);

                    simdgeo::geom_h_y(buffer, 1550, 752, 860, 1, 3, ncols, beta);

                    simdgeo::geom_h_z(buffer, 1613, 752, 860, 1, 3, ncols, beta);

                    simdgeo::geom_i_x(buffer, 1676, 797, 944, 1, 3, ncols, beta);

                    simdgeo::geom_i_y(buffer, 1760, 797, 944, 1, 3, ncols, beta);

                    simdgeo::geom_i_z(buffer, 1844, 797, 944, 1, 3, ncols, beta);

                    simdgeo::geom_k_x(buffer, 1928, 860, 1052, 1, 3, ncols, beta);

                    simdgeo::geom_k_y(buffer, 2036, 860, 1052, 1, 3, ncols, beta);

                    simdgeo::geom_k_z(buffer, 2144, 860, 1052, 1, 3, ncols, beta);

                    simdgeo::geom_l_x(buffer, 2252, 944, 1187, 1, 3, ncols, beta);

                    simdgeo::geom_l_y(buffer, 2387, 944, 1187, 1, 3, ncols, beta);

                    simdgeo::geom_l_z(buffer, 2522, 944, 1187, 1, 3, ncols, beta);

                    simdfunc::contract_primitives(buffer, 2657, 1352, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2747, 1397, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2837, 1442, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2927, 752, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3017, 1487, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3143, 1550, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3269, 1613, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3395, 797, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3521, 1676, 84, ncols);

                    simdfunc::contract_primitives(buffer, 3689, 1760, 84, ncols);

                    simdfunc::contract_primitives(buffer, 3857, 1844, 84, ncols);

                    simdfunc::contract_primitives(buffer, 4025, 860, 84, ncols);

                    simdfunc::contract_primitives(buffer, 4193, 1928, 108, ncols);

                    simdfunc::contract_primitives(buffer, 4409, 2036, 108, ncols);

                    simdfunc::contract_primitives(buffer, 4625, 2144, 108, ncols);

                    simdfunc::contract_primitives(buffer, 4841, 944, 108, ncols);

                    simdfunc::contract_primitives(buffer, 5057, 2252, 135, ncols);

                    simdfunc::contract_primitives(buffer, 5327, 2387, 135, ncols);

                    simdfunc::contract_primitives(buffer, 5597, 2522, 135, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 2702, 2657, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2792, 2747, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2882, 2837, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2972, 2927, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3080, 3017, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3206, 3143, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3332, 3269, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3458, 3395, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3605, 3521, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3773, 3689, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3941, 3857, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4109, 4025, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4301, 4193, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4517, 4409, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4733, 4625, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4949, 4841, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5192, 5057, 45, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5462, 5327, 45, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5732, 5597, 45, 1, nmax);

        simdtrf::compute_hrr_geom_010x_pg(buffer, coordinates, 5867, 2702, 2972, 3080, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pg(buffer, coordinates, 6002, 2792, 2972, 3206, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pg(buffer, coordinates, 6137, 2882, 2972, 3332, 3, nmax);

        simdtrf::compute_hrr_pg(buffer, coordinates, 6272, 2972, 3458, 3, nmax);

        simdtrf::compute_hrr_geom_010x_ph(buffer, coordinates, 6407, 3080, 3458, 3605, 3, nmax);

        simdtrf::compute_hrr_geom_010y_ph(buffer, coordinates, 6596, 3206, 3458, 3773, 3, nmax);

        simdtrf::compute_hrr_geom_010z_ph(buffer, coordinates, 6785, 3332, 3458, 3941, 3, nmax);

        simdtrf::compute_hrr_ph(buffer, coordinates, 6974, 3458, 4109, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pi(buffer, coordinates, 7163, 3605, 4109, 4301, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pi(buffer, coordinates, 7415, 3773, 4109, 4517, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pi(buffer, coordinates, 7667, 3941, 4109, 4733, 3, nmax);

        simdtrf::compute_hrr_pi(buffer, coordinates, 7919, 4109, 4949, 3, nmax);

        simdtrf::compute_hrr_geom_010x_pk(buffer, coordinates, 8171, 4301, 4949, 5192, 3, nmax);

        simdtrf::compute_hrr_geom_010y_pk(buffer, coordinates, 8495, 4517, 4949, 5462, 3, nmax);

        simdtrf::compute_hrr_geom_010z_pk(buffer, coordinates, 8819, 4733, 4949, 5732, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dg(buffer, coordinates, 9143, 5867, 6272, 6407, 3, nmax);

        simdtrf::compute_hrr_geom_010y_dg(buffer, coordinates, 9413, 6002, 6272, 6596, 3, nmax);

        simdtrf::compute_hrr_geom_010z_dg(buffer, coordinates, 9683, 6137, 6272, 6785, 3, nmax);

        simdtrf::compute_hrr_dg(buffer, coordinates, 9953, 6272, 6974, 3, nmax);

        simdtrf::compute_hrr_geom_010x_dh(buffer, coordinates, 10223, 6407, 6974, 7163, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_dh(buffer, coordinates, 10601, 6596, 6974, 7415, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_dh(buffer, coordinates, 10979, 6785, 6974, 7667, 3,
                                          nmax);

        simdtrf::compute_hrr_dh(buffer, coordinates, 11357, 6974, 7919, 3, nmax);

        simdtrf::compute_hrr_geom_010x_di(buffer, coordinates, 11735, 7163, 7919, 8171, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_di(buffer, coordinates, 12239, 7415, 7919, 8495, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_di(buffer, coordinates, 12743, 7667, 7919, 8819, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_fg(buffer, coordinates, 13247, 9143, 9953, 10223, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fg(buffer, coordinates, 13697, 9413, 9953, 10601, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fg(buffer, coordinates, 14147, 9683, 9953, 10979, 3,
                                          nmax);

        simdtrf::compute_hrr_fg(buffer, coordinates, 14597, 9953, 11357, 3, nmax);

        simdtrf::compute_hrr_geom_010x_fh(buffer, coordinates, 15047, 10223, 11357, 11735, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_fh(buffer, coordinates, 15677, 10601, 11357, 12239, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_fh(buffer, coordinates, 16307, 10979, 11357, 12743, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010x_gg(buffer, coordinates, 16937, 13247, 14597, 15047, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010y_gg(buffer, coordinates, 17612, 13697, 14597, 15677, 3,
                                          nmax);

        simdtrf::compute_hrr_geom_010z_gg(buffer, coordinates, 18287, 14147, 14597, 16307, 3,
                                          nmax);

        simdtrf::transform_g_inner(buffer, 18962, 16937, 15, 3, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 18962, 27, nmax);

        simdtrf::transform_g_inner(buffer, 18962, 17612, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 243 * nvalues + n * npairs, nvalues, buffer, 18962,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 18962, 18287, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 486 * nvalues + n * npairs, nvalues, buffer, 18962,
                                   27, nmax);
    }

    for (size_t m = 0; m < 729; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
