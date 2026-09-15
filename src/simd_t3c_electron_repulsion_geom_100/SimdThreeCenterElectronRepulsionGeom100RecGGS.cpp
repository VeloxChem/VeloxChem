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


#include "SimdThreeCenterElectronRepulsionGeom100RecGGS.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecMSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdTransferGD.hpp"
#include "SimdTransferGF.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferGeom100XGD.hpp"
#include "SimdTransferGeom100XGF.hpp"
#include "SimdTransferGeom100XGG.hpp"
#include "SimdTransferGeom100XGP.hpp"
#include "SimdTransferGeom100XHD.hpp"
#include "SimdTransferGeom100XHF.hpp"
#include "SimdTransferGeom100XHP.hpp"
#include "SimdTransferGeom100XID.hpp"
#include "SimdTransferGeom100XIP.hpp"
#include "SimdTransferGeom100XKP.hpp"
#include "SimdTransferGeom100YGD.hpp"
#include "SimdTransferGeom100YGF.hpp"
#include "SimdTransferGeom100YGG.hpp"
#include "SimdTransferGeom100YGP.hpp"
#include "SimdTransferGeom100YHD.hpp"
#include "SimdTransferGeom100YHF.hpp"
#include "SimdTransferGeom100YHP.hpp"
#include "SimdTransferGeom100YID.hpp"
#include "SimdTransferGeom100YIP.hpp"
#include "SimdTransferGeom100YKP.hpp"
#include "SimdTransferGeom100ZGD.hpp"
#include "SimdTransferGeom100ZGF.hpp"
#include "SimdTransferGeom100ZGG.hpp"
#include "SimdTransferGeom100ZGP.hpp"
#include "SimdTransferGeom100ZHD.hpp"
#include "SimdTransferGeom100ZHF.hpp"
#include "SimdTransferGeom100ZHP.hpp"
#include "SimdTransferGeom100ZID.hpp"
#include "SimdTransferGeom100ZIP.hpp"
#include "SimdTransferGeom100ZKP.hpp"
#include "SimdTransferHD.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransferIP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_ggs_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_ggs_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 6727, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 243 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 6727, 1157, 1025, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 9, ncols,
                                                             fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 17, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 20, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 23, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 7, 8,
                                                                       17, 20, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 8, 9,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 9, 10,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 10, 11,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 11, 12,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 12, 13,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 13, 14,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 14, 15,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 17, 20,
                                                                       44, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 20, 23,
                                                                       50, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 23, 26,
                                                                       56, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 26, 29,
                                                                       62, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 29, 32,
                                                                       68, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 32, 35,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 35, 38,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 162, 0, 3, 44, 50,
                                                                       92, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 177, 0, 3, 50, 56,
                                                                       102, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 56, 62,
                                                                       112, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 207, 0, 3, 62, 68,
                                                                       122, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 222, 0, 3, 68, 74,
                                                                       132, 142, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 237, 0, 3, 74, 80,
                                                                       142, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 252, 0, 3, 92,
                                                                       102, 162, 177, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 102,
                                                                       112, 177, 192, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 294, 0, 3, 112,
                                                                       122, 192, 207, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 315, 0, 3, 122,
                                                                       132, 207, 222, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 336, 0, 3, 132,
                                                                       142, 222, 237, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 357, 0, 3, 162,
                                                                       177, 252, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 385, 0, 3, 177,
                                                                       192, 273, 294, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 413, 0, 3, 192,
                                                                       207, 294, 315, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 441, 0, 3, 207,
                                                                       222, 315, 336, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 469, 0, 3, 252,
                                                                       273, 357, 385, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 505, 0, 3, 273,
                                                                       294, 385, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 541, 0, 3, 294,
                                                                       315, 413, 441, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 577, 0, 3, 357,
                                                                       385, 469, 505, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 622, 0, 3, 385,
                                                                       413, 505, 541, ncols,
                                                                       gamma, p, q);

                    compute_prim_mss_three_center_electron_repulsion_0(buffer, 667, 0, 3, 469,
                                                                       505, 577, 622, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_g_x(buffer, 722, 92, 252, 1, 1, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 737, 92, 252, 1, 1, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 752, 92, 252, 1, 1, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 767, 162, 357, 1, 1, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 788, 162, 357, 1, 1, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 809, 162, 357, 1, 1, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 830, 252, 469, 1, 1, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 858, 252, 469, 1, 1, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 886, 252, 469, 1, 1, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 914, 357, 577, 1, 1, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 950, 357, 577, 1, 1, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 986, 357, 577, 1, 1, ncols, alpha);

                    simdgeo::geom_l_x(buffer, 1022, 469, 667, 1, 1, ncols, alpha);

                    simdgeo::geom_l_y(buffer, 1067, 469, 667, 1, 1, ncols, alpha);

                    simdgeo::geom_l_z(buffer, 1112, 469, 667, 1, 1, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 1157, 722, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1187, 737, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1217, 752, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1247, 162, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1277, 767, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1319, 788, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1361, 809, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1403, 252, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1445, 830, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1501, 858, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1557, 886, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1613, 357, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1669, 914, 36, ncols);

                    simdfunc::contract_primitives(buffer, 1741, 950, 36, ncols);

                    simdfunc::contract_primitives(buffer, 1813, 986, 36, ncols);

                    simdfunc::contract_primitives(buffer, 1885, 469, 36, ncols);

                    simdfunc::contract_primitives(buffer, 1957, 1022, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2047, 1067, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2137, 1112, 45, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 1172, 1157, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1202, 1187, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1232, 1217, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1262, 1247, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1298, 1277, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1340, 1319, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1382, 1361, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1424, 1403, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1473, 1445, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1529, 1501, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1585, 1557, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1641, 1613, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1705, 1669, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1777, 1741, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1849, 1813, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1921, 1885, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2002, 1957, 45, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2092, 2047, 45, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2182, 2137, 45, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 2227, 1172, 1262,
                                                       1298, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 2272, 1202, 1262,
                                                       1340, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 2317, 1232, 1262,
                                                       1382, 1, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 2362, 1262, 1424, 1, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 2407, 1298, 1424,
                                                       1473, 1, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 2470, 1340, 1424,
                                                       1529, 1, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 2533, 1382, 1424,
                                                       1585, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 2596, 1424, 1641, 1, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 2659, 1473, 1641,
                                                       1705, 1, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 2743, 1529, 1641,
                                                       1777, 1, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 2827, 1585, 1641,
                                                       1849, 1, nmax);

        simdtrf::compute_hrr_ip_out_of_first(buffer, coordinates, 2911, 1641, 1921, 1, nmax);

        simdtrf::compute_hrr_geom_100x_kp_out_of_first(buffer, coordinates, 2995, 1705, 1921,
                                                       2002, 1, nmax);

        simdtrf::compute_hrr_geom_100y_kp_out_of_first(buffer, coordinates, 3103, 1777, 1921,
                                                       2092, 1, nmax);

        simdtrf::compute_hrr_geom_100z_kp_out_of_first(buffer, coordinates, 3211, 1849, 1921,
                                                       2182, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 3319, 2227, 2362,
                                                       2407, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 3409, 2272, 2362,
                                                       2470, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 3499, 2317, 2362,
                                                       2533, 1, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 3589, 2362, 2596, 1, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 3679, 2407, 2596,
                                                       2659, 1, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 3805, 2470, 2596,
                                                       2743, 1, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 3931, 2533, 2596,
                                                       2827, 1, nmax);

        simdtrf::compute_hrr_hd_out_of_first(buffer, coordinates, 4057, 2596, 2911, 1, nmax);

        simdtrf::compute_hrr_geom_100x_id_out_of_first(buffer, coordinates, 4183, 2659, 2911,
                                                       2995, 1, nmax);

        simdtrf::compute_hrr_geom_100y_id_out_of_first(buffer, coordinates, 4351, 2743, 2911,
                                                       3103, 1, nmax);

        simdtrf::compute_hrr_geom_100z_id_out_of_first(buffer, coordinates, 4519, 2827, 2911,
                                                       3211, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 4687, 3319, 3589,
                                                       3679, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 4837, 3409, 3589,
                                                       3805, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 4987, 3499, 3589,
                                                       3931, 1, nmax);

        simdtrf::compute_hrr_gf_out_of_first(buffer, coordinates, 5137, 3589, 4057, 1, nmax);

        simdtrf::compute_hrr_geom_100x_hf_out_of_first(buffer, coordinates, 5287, 3679, 4057,
                                                       4183, 1, nmax);

        simdtrf::compute_hrr_geom_100y_hf_out_of_first(buffer, coordinates, 5497, 3805, 4057,
                                                       4351, 1, nmax);

        simdtrf::compute_hrr_geom_100z_hf_out_of_first(buffer, coordinates, 5707, 3931, 4057,
                                                       4519, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gg_out_of_first(buffer, coordinates, 5917, 4687, 5137,
                                                       5287, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gg_out_of_first(buffer, coordinates, 6142, 4837, 5137,
                                                       5497, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gg_out_of_first(buffer, coordinates, 6367, 4987, 5137,
                                                       5707, 1, nmax);

        simdtrf::transform_g_inner(buffer, 6592, 5917, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 6592, 9, nmax);

        simdtrf::transform_g_inner(buffer, 6592, 6142, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 81 * nvalues + n * npairs, nvalues, buffer, 6592, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 6592, 6367, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 162 * nvalues + n * npairs, nvalues, buffer, 6592, 9,
                                   nmax);
    }

    for (size_t m = 0; m < 243; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
