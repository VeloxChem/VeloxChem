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


#include "SimdThreeCenterElectronRepulsionGeom100RecPGP.hpp"

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
#include "SimdGeometryP1.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDD.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGeom100XDD.hpp"
#include "SimdTransferGeom100XDF.hpp"
#include "SimdTransferGeom100XDP.hpp"
#include "SimdTransferGeom100XFD.hpp"
#include "SimdTransferGeom100XFP.hpp"
#include "SimdTransferGeom100XGP.hpp"
#include "SimdTransferGeom100XPD.hpp"
#include "SimdTransferGeom100XPF.hpp"
#include "SimdTransferGeom100XPG.hpp"
#include "SimdTransferGeom100XPP.hpp"
#include "SimdTransferGeom100YDD.hpp"
#include "SimdTransferGeom100YDF.hpp"
#include "SimdTransferGeom100YDP.hpp"
#include "SimdTransferGeom100YFD.hpp"
#include "SimdTransferGeom100YFP.hpp"
#include "SimdTransferGeom100YGP.hpp"
#include "SimdTransferGeom100YPD.hpp"
#include "SimdTransferGeom100YPF.hpp"
#include "SimdTransferGeom100YPG.hpp"
#include "SimdTransferGeom100YPP.hpp"
#include "SimdTransferGeom100ZDD.hpp"
#include "SimdTransferGeom100ZDF.hpp"
#include "SimdTransferGeom100ZDP.hpp"
#include "SimdTransferGeom100ZFD.hpp"
#include "SimdTransferGeom100ZFP.hpp"
#include "SimdTransferGeom100ZGP.hpp"
#include "SimdTransferGeom100ZPD.hpp"
#include "SimdTransferGeom100ZPF.hpp"
#include "SimdTransferGeom100ZPG.hpp"
#include "SimdTransferGeom100ZPP.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransferPF.hpp"
#include "SimdTransferPP.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_pgp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_pgp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 5821, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 5821, 964, 1131, dimensions);

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
                                                        5, 6, 7}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 14, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 17, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 20, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 23, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 7, 8,
                                                                       14, 17, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 8, 9,
                                                                       17, 20, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 9, 10,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 10, 11,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 11, 12,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 14, 17,
                                                                       32, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 17, 20,
                                                                       38, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 20, 23,
                                                                       44, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 23, 26,
                                                                       50, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 32, 38,
                                                                       62, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 117, 0, 3, 38, 44,
                                                                       72, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 44, 50,
                                                                       82, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 147, 0, 3, 62, 72,
                                                                       102, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 168, 0, 3, 72, 82,
                                                                       117, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 189, 0, 3, 102,
                                                                       117, 147, 168, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 217, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 220, 3, 7, 14,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 229, 3, 14, 32,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 247, 3, 32, 62,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 277, 3, 62, 102,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 322, 3, 102, 147,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 385, 3, 147, 189,
                                                                       ncols, p, q);

                    simdgeo::geom_p_x(buffer, 469, 217, 229, 1, 3, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 478, 217, 229, 1, 3, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 487, 217, 229, 1, 3, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 496, 220, 247, 1, 3, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 514, 220, 247, 1, 3, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 532, 220, 247, 1, 3, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 550, 229, 277, 1, 3, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 580, 229, 277, 1, 3, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 610, 229, 277, 1, 3, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 640, 247, 322, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 685, 247, 322, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 730, 247, 322, 1, 3, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 775, 277, 385, 1, 3, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 838, 277, 385, 1, 3, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 901, 277, 385, 1, 3, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 964, 469, 9, ncols);

                    simdfunc::contract_primitives(buffer, 982, 478, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1000, 487, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1018, 220, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1036, 496, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1072, 514, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1108, 532, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1144, 229, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1180, 550, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1240, 580, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1300, 610, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1360, 247, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1420, 640, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1510, 685, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1600, 730, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1690, 277, 45, ncols);

                    simdfunc::contract_primitives(buffer, 1780, 775, 63, ncols);

                    simdfunc::contract_primitives(buffer, 1906, 838, 63, ncols);

                    simdfunc::contract_primitives(buffer, 2032, 901, 63, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 973, 964, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 991, 982, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1009, 1000, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1027, 1018, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1054, 1036, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1090, 1072, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1126, 1108, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1162, 1144, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1210, 1180, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1270, 1240, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1330, 1300, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1390, 1360, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1465, 1420, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1555, 1510, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1645, 1600, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1735, 1690, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1843, 1780, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1969, 1906, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2095, 2032, 21, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 2158, 973, 1027,
                                                       1054, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 2185, 991, 1027,
                                                       1090, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 2212, 1009, 1027,
                                                       1126, 3, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 2239, 1027, 1162, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 2266, 1054, 1162,
                                                       1210, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 2320, 1090, 1162,
                                                       1270, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 2374, 1126, 1162,
                                                       1330, 3, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 2428, 1162, 1390, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 2482, 1210, 1390,
                                                       1465, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 2572, 1270, 1390,
                                                       1555, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 2662, 1330, 1390,
                                                       1645, 3, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 2752, 1390, 1735, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 2842, 1465, 1735,
                                                       1843, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 2977, 1555, 1735,
                                                       1969, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 3112, 1645, 1735,
                                                       2095, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 3247, 2158, 2239,
                                                       2266, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 3301, 2185, 2239,
                                                       2320, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 3355, 2212, 2239,
                                                       2374, 3, nmax);

        simdtrf::compute_hrr_pd_out_of_first(buffer, coordinates, 3409, 2239, 2428, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 3463, 2266, 2428,
                                                       2482, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 3571, 2320, 2428,
                                                       2572, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 3679, 2374, 2428,
                                                       2662, 3, nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 3787, 2428, 2752, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 3895, 2482, 2752,
                                                       2842, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 4075, 2572, 2752,
                                                       2977, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 4255, 2662, 2752,
                                                       3112, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pf_out_of_first(buffer, coordinates, 4435, 3247, 3409,
                                                       3463, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pf_out_of_first(buffer, coordinates, 4525, 3301, 3409,
                                                       3571, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pf_out_of_first(buffer, coordinates, 4615, 3355, 3409,
                                                       3679, 3, nmax);

        simdtrf::compute_hrr_pf_out_of_first(buffer, coordinates, 4705, 3409, 3787, 3, nmax);

        simdtrf::compute_hrr_geom_100x_df_out_of_first(buffer, coordinates, 4795, 3463, 3787,
                                                       3895, 3, nmax);

        simdtrf::compute_hrr_geom_100y_df_out_of_first(buffer, coordinates, 4975, 3571, 3787,
                                                       4075, 3, nmax);

        simdtrf::compute_hrr_geom_100z_df_out_of_first(buffer, coordinates, 5155, 3679, 3787,
                                                       4255, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pg_out_of_first(buffer, coordinates, 5335, 4435, 4705,
                                                       4795, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pg_out_of_first(buffer, coordinates, 5470, 4525, 4705,
                                                       4975, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pg_out_of_first(buffer, coordinates, 5605, 4615, 4705,
                                                       5155, 3, nmax);

        simdtrf::transform_g_inner(buffer, 5740, 5335, 3, 3, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 5740, 27, nmax);

        simdtrf::transform_g_inner(buffer, 5740, 5470, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 81 * nvalues + n * npairs, nvalues, buffer, 5740, 27,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 5740, 5605, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 162 * nvalues + n * npairs, nvalues, buffer, 5740,
                                   27, nmax);
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
