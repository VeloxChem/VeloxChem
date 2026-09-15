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


#include "SimdThreeCenterElectronRepulsionGeom100RecDGP.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdTransferDD.hpp"
#include "SimdTransferDF.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferFD.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferGeom100XDD.hpp"
#include "SimdTransferGeom100XDF.hpp"
#include "SimdTransferGeom100XDG.hpp"
#include "SimdTransferGeom100XDP.hpp"
#include "SimdTransferGeom100XFD.hpp"
#include "SimdTransferGeom100XFF.hpp"
#include "SimdTransferGeom100XFP.hpp"
#include "SimdTransferGeom100XGD.hpp"
#include "SimdTransferGeom100XGP.hpp"
#include "SimdTransferGeom100XHP.hpp"
#include "SimdTransferGeom100YDD.hpp"
#include "SimdTransferGeom100YDF.hpp"
#include "SimdTransferGeom100YDG.hpp"
#include "SimdTransferGeom100YDP.hpp"
#include "SimdTransferGeom100YFD.hpp"
#include "SimdTransferGeom100YFF.hpp"
#include "SimdTransferGeom100YFP.hpp"
#include "SimdTransferGeom100YGD.hpp"
#include "SimdTransferGeom100YGP.hpp"
#include "SimdTransferGeom100YHP.hpp"
#include "SimdTransferGeom100ZDD.hpp"
#include "SimdTransferGeom100ZDF.hpp"
#include "SimdTransferGeom100ZDG.hpp"
#include "SimdTransferGeom100ZDP.hpp"
#include "SimdTransferGeom100ZFD.hpp"
#include "SimdTransferGeom100ZFF.hpp"
#include "SimdTransferGeom100ZFP.hpp"
#include "SimdTransferGeom100ZGD.hpp"
#include "SimdTransferGeom100ZGP.hpp"
#include "SimdTransferGeom100ZHP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_dgp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_dgp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 9403, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 405 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 9403, 1414, 1668, dimensions);

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 6, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 15, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 18, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 21, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 24, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 27, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 7, 8,
                                                                       15, 18, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 8, 9,
                                                                       18, 21, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 9, 10,
                                                                       21, 24, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 10, 11,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 60, 0, 3, 11, 12,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 66, 0, 3, 12, 13,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 15, 18,
                                                                       36, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 18, 21,
                                                                       42, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 21, 24,
                                                                       48, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 24, 27,
                                                                       54, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 27, 30,
                                                                       60, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 36, 42,
                                                                       72, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 137, 0, 3, 42, 48,
                                                                       82, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 48, 54,
                                                                       92, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 167, 0, 3, 54, 60,
                                                                       102, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 182, 0, 3, 72, 82,
                                                                       122, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 82, 92,
                                                                       137, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 224, 0, 3, 92,
                                                                       102, 152, 167, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 245, 0, 3, 122,
                                                                       137, 182, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 137,
                                                                       152, 203, 224, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 301, 0, 3, 182,
                                                                       203, 245, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 337, 3, 7, 15,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 346, 3, 15, 36,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 364, 3, 36, 72,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 394, 3, 72, 122,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 439, 3, 122, 182,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 502, 3, 182, 245,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 586, 3, 245, 301,
                                                                       ncols, p, q);

                    simdgeo::geom_d_x(buffer, 694, 337, 364, 1, 3, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 712, 337, 364, 1, 3, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 730, 337, 364, 1, 3, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 748, 346, 394, 1, 3, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 778, 346, 394, 1, 3, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 808, 346, 394, 1, 3, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 838, 364, 439, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 883, 364, 439, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 928, 364, 439, 1, 3, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 973, 394, 502, 1, 3, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 1036, 394, 502, 1, 3, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 1099, 394, 502, 1, 3, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 1162, 439, 586, 1, 3, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 1246, 439, 586, 1, 3, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 1330, 439, 586, 1, 3, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 1414, 694, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1450, 712, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1486, 730, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1522, 346, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1558, 748, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1618, 778, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1678, 808, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1738, 364, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1798, 838, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1888, 883, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1978, 928, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2068, 394, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2158, 973, 63, ncols);

                    simdfunc::contract_primitives(buffer, 2284, 1036, 63, ncols);

                    simdfunc::contract_primitives(buffer, 2410, 1099, 63, ncols);

                    simdfunc::contract_primitives(buffer, 2536, 439, 63, ncols);

                    simdfunc::contract_primitives(buffer, 2662, 1162, 84, ncols);

                    simdfunc::contract_primitives(buffer, 2830, 1246, 84, ncols);

                    simdfunc::contract_primitives(buffer, 2998, 1330, 84, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 1432, 1414, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1468, 1450, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1504, 1486, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1540, 1522, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1588, 1558, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1648, 1618, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1708, 1678, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1768, 1738, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1843, 1798, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1933, 1888, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2023, 1978, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2113, 2068, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2221, 2158, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2347, 2284, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2473, 2410, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2599, 2536, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2746, 2662, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2914, 2830, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3082, 2998, 28, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 3166, 1432, 1540,
                                                       1588, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 3220, 1468, 1540,
                                                       1648, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 3274, 1504, 1540,
                                                       1708, 3, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 3328, 1540, 1768, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 3382, 1588, 1768,
                                                       1843, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 3472, 1648, 1768,
                                                       1933, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 3562, 1708, 1768,
                                                       2023, 3, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 3652, 1768, 2113, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 3742, 1843, 2113,
                                                       2221, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 3877, 1933, 2113,
                                                       2347, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 4012, 2023, 2113,
                                                       2473, 3, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 4147, 2113, 2599, 3, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 4282, 2221, 2599,
                                                       2746, 3, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 4471, 2347, 2599,
                                                       2914, 3, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 4660, 2473, 2599,
                                                       3082, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 4849, 3166, 3328,
                                                       3382, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 4957, 3220, 3328,
                                                       3472, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 5065, 3274, 3328,
                                                       3562, 3, nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 5173, 3328, 3652, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 5281, 3382, 3652,
                                                       3742, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 5461, 3472, 3652,
                                                       3877, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 5641, 3562, 3652,
                                                       4012, 3, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 5821, 3652, 4147, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 6001, 3742, 4147,
                                                       4282, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 6271, 3877, 4147,
                                                       4471, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 6541, 4012, 4147,
                                                       4660, 3, nmax);

        simdtrf::compute_hrr_geom_100x_df_out_of_first(buffer, coordinates, 6811, 4849, 5173,
                                                       5281, 3, nmax);

        simdtrf::compute_hrr_geom_100y_df_out_of_first(buffer, coordinates, 6991, 4957, 5173,
                                                       5461, 3, nmax);

        simdtrf::compute_hrr_geom_100z_df_out_of_first(buffer, coordinates, 7171, 5065, 5173,
                                                       5641, 3, nmax);

        simdtrf::compute_hrr_df_out_of_first(buffer, coordinates, 7351, 5173, 5821, 3, nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 7531, 5281, 5821,
                                                       6001, 3, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 7831, 5461, 5821,
                                                       6271, 3, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 8131, 5641, 5821,
                                                       6541, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dg_out_of_first(buffer, coordinates, 8431, 6811, 7351,
                                                       7531, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dg_out_of_first(buffer, coordinates, 8701, 6991, 7351,
                                                       7831, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dg_out_of_first(buffer, coordinates, 8971, 7171, 7351,
                                                       8131, 3, nmax);

        simdtrf::transform_g_inner(buffer, 9241, 8431, 6, 3, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 9241, 27, nmax);

        simdtrf::transform_g_inner(buffer, 9241, 8701, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 135 * nvalues + n * npairs, nvalues, buffer, 9241,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 9241, 8971, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 270 * nvalues + n * npairs, nvalues, buffer, 9241,
                                   27, nmax);
    }

    for (size_t m = 0; m < 405; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
