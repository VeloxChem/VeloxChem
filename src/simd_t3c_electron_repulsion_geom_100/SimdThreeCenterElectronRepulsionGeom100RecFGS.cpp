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


#include "SimdThreeCenterElectronRepulsionGeom100RecFGS.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdTransferFD.hpp"
#include "SimdTransferFF.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGD.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferGeom100XFD.hpp"
#include "SimdTransferGeom100XFF.hpp"
#include "SimdTransferGeom100XFG.hpp"
#include "SimdTransferGeom100XFP.hpp"
#include "SimdTransferGeom100XGD.hpp"
#include "SimdTransferGeom100XGF.hpp"
#include "SimdTransferGeom100XGP.hpp"
#include "SimdTransferGeom100XHD.hpp"
#include "SimdTransferGeom100XHP.hpp"
#include "SimdTransferGeom100XIP.hpp"
#include "SimdTransferGeom100YFD.hpp"
#include "SimdTransferGeom100YFF.hpp"
#include "SimdTransferGeom100YFG.hpp"
#include "SimdTransferGeom100YFP.hpp"
#include "SimdTransferGeom100YGD.hpp"
#include "SimdTransferGeom100YGF.hpp"
#include "SimdTransferGeom100YGP.hpp"
#include "SimdTransferGeom100YHD.hpp"
#include "SimdTransferGeom100YHP.hpp"
#include "SimdTransferGeom100YIP.hpp"
#include "SimdTransferGeom100ZFD.hpp"
#include "SimdTransferGeom100ZFF.hpp"
#include "SimdTransferGeom100ZFG.hpp"
#include "SimdTransferGeom100ZFP.hpp"
#include "SimdTransferGeom100ZGD.hpp"
#include "SimdTransferGeom100ZGF.hpp"
#include "SimdTransferGeom100ZGP.hpp"
#include "SimdTransferGeom100ZHD.hpp"
#include "SimdTransferGeom100ZHP.hpp"
#include "SimdTransferGeom100ZIP.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_fgs_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_fgs_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 4812, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 4812, 832, 772, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 8, ncols,
                                                             fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 16, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 19, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 22, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 7, 8,
                                                                       16, 19, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 8, 9,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 9, 10,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 10, 11,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 11, 12,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 12, 13,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 76, 0, 3, 13, 14,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 16, 19,
                                                                       40, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 19, 22,
                                                                       46, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 22, 25,
                                                                       52, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 25, 28,
                                                                       58, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 28, 31,
                                                                       64, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 31, 34,
                                                                       70, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 142, 0, 3, 40, 46,
                                                                       82, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 157, 0, 3, 46, 52,
                                                                       92, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 172, 0, 3, 52, 58,
                                                                       102, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 187, 0, 3, 58, 64,
                                                                       112, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 202, 0, 3, 64, 70,
                                                                       122, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 217, 0, 3, 82, 92,
                                                                       142, 157, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 92,
                                                                       102, 157, 172, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 259, 0, 3, 102,
                                                                       112, 172, 187, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 280, 0, 3, 112,
                                                                       122, 187, 202, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 301, 0, 3, 142,
                                                                       157, 217, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 329, 0, 3, 157,
                                                                       172, 238, 259, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 357, 0, 3, 172,
                                                                       187, 259, 280, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 385, 0, 3, 217,
                                                                       238, 301, 329, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 421, 0, 3, 238,
                                                                       259, 329, 357, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 457, 0, 3, 301,
                                                                       329, 385, 421, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_f_x(buffer, 502, 40, 142, 1, 1, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 512, 40, 142, 1, 1, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 522, 40, 142, 1, 1, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 532, 82, 217, 1, 1, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 547, 82, 217, 1, 1, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 562, 82, 217, 1, 1, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 577, 142, 301, 1, 1, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 598, 142, 301, 1, 1, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 619, 142, 301, 1, 1, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 640, 217, 385, 1, 1, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 668, 217, 385, 1, 1, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 696, 217, 385, 1, 1, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 724, 301, 457, 1, 1, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 760, 301, 457, 1, 1, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 796, 301, 457, 1, 1, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 832, 502, 10, ncols);

                    simdfunc::contract_primitives(buffer, 852, 512, 10, ncols);

                    simdfunc::contract_primitives(buffer, 872, 522, 10, ncols);

                    simdfunc::contract_primitives(buffer, 892, 82, 10, ncols);

                    simdfunc::contract_primitives(buffer, 912, 532, 15, ncols);

                    simdfunc::contract_primitives(buffer, 942, 547, 15, ncols);

                    simdfunc::contract_primitives(buffer, 972, 562, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1002, 142, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1032, 577, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1074, 598, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1116, 619, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1158, 217, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1200, 640, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1256, 668, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1312, 696, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1368, 301, 28, ncols);

                    simdfunc::contract_primitives(buffer, 1424, 724, 36, ncols);

                    simdfunc::contract_primitives(buffer, 1496, 760, 36, ncols);

                    simdfunc::contract_primitives(buffer, 1568, 796, 36, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 842, 832, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 862, 852, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 882, 872, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 902, 892, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 927, 912, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 957, 942, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 987, 972, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1017, 1002, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1053, 1032, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1095, 1074, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1137, 1116, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1179, 1158, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1228, 1200, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1284, 1256, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1340, 1312, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1396, 1368, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1460, 1424, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1532, 1496, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1604, 1568, 36, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 1640, 842, 902, 927,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 1670, 862, 902, 957,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 1700, 882, 902, 987,
                                                       1, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 1730, 902, 1017, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 1760, 927, 1017,
                                                       1053, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 1805, 957, 1017,
                                                       1095, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 1850, 987, 1017,
                                                       1137, 1, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 1895, 1017, 1179, 1, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 1940, 1053, 1179,
                                                       1228, 1, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 2003, 1095, 1179,
                                                       1284, 1, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 2066, 1137, 1179,
                                                       1340, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 2129, 1179, 1396, 1, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 2192, 1228, 1396,
                                                       1460, 1, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 2276, 1284, 1396,
                                                       1532, 1, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 2360, 1340, 1396,
                                                       1604, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 2444, 1640, 1730,
                                                       1760, 1, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 2504, 1670, 1730,
                                                       1805, 1, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 2564, 1700, 1730,
                                                       1850, 1, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 2624, 1730, 1895, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 2684, 1760, 1895,
                                                       1940, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 2774, 1805, 1895,
                                                       2003, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 2864, 1850, 1895,
                                                       2066, 1, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 2954, 1895, 2129, 1, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 3044, 1940, 2129,
                                                       2192, 1, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 3170, 2003, 2129,
                                                       2276, 1, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 3296, 2066, 2129,
                                                       2360, 1, nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 3422, 2444, 2624,
                                                       2684, 1, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 3522, 2504, 2624,
                                                       2774, 1, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 3622, 2564, 2624,
                                                       2864, 1, nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 3722, 2624, 2954, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 3822, 2684, 2954,
                                                       3044, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 3972, 2774, 2954,
                                                       3170, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 4122, 2864, 2954,
                                                       3296, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fg_out_of_first(buffer, coordinates, 4272, 3422, 3722,
                                                       3822, 1, nmax);

        simdtrf::compute_hrr_geom_100y_fg_out_of_first(buffer, coordinates, 4422, 3522, 3722,
                                                       3972, 1, nmax);

        simdtrf::compute_hrr_geom_100z_fg_out_of_first(buffer, coordinates, 4572, 3622, 3722,
                                                       4122, 1, nmax);

        simdtrf::transform_g_inner(buffer, 4722, 4272, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 4722, 9, nmax);

        simdtrf::transform_g_inner(buffer, 4722, 4422, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 4722, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 4722, 4572, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 126 * nvalues + n * npairs, nvalues, buffer, 4722, 9,
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
