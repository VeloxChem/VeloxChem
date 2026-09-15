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


#include "SimdThreeCenterElectronRepulsionGeom100RecPGS.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
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
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_geom_100_pgs_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_100_pgs_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2001, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 81 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 2001, 382, 377, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 6, 3, 6, ncols,
                                                             fj, i * nprim_b + j, fq);

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

                    simdgeo::geom_p_x(buffer, 217, 7, 32, 1, 1, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 220, 7, 32, 1, 1, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 223, 7, 32, 1, 1, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 226, 14, 62, 1, 1, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 232, 14, 62, 1, 1, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 238, 14, 62, 1, 1, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 244, 32, 102, 1, 1, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 254, 32, 102, 1, 1, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 264, 32, 102, 1, 1, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 274, 62, 147, 1, 1, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 289, 62, 147, 1, 1, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 304, 62, 147, 1, 1, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 319, 102, 189, 1, 1, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 340, 102, 189, 1, 1, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 361, 102, 189, 1, 1, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 382, 217, 3, ncols);

                    simdfunc::contract_primitives(buffer, 388, 220, 3, ncols);

                    simdfunc::contract_primitives(buffer, 394, 223, 3, ncols);

                    simdfunc::contract_primitives(buffer, 400, 14, 3, ncols);

                    simdfunc::contract_primitives(buffer, 406, 226, 6, ncols);

                    simdfunc::contract_primitives(buffer, 418, 232, 6, ncols);

                    simdfunc::contract_primitives(buffer, 430, 238, 6, ncols);

                    simdfunc::contract_primitives(buffer, 442, 32, 6, ncols);

                    simdfunc::contract_primitives(buffer, 454, 244, 10, ncols);

                    simdfunc::contract_primitives(buffer, 474, 254, 10, ncols);

                    simdfunc::contract_primitives(buffer, 494, 264, 10, ncols);

                    simdfunc::contract_primitives(buffer, 514, 62, 10, ncols);

                    simdfunc::contract_primitives(buffer, 534, 274, 15, ncols);

                    simdfunc::contract_primitives(buffer, 564, 289, 15, ncols);

                    simdfunc::contract_primitives(buffer, 594, 304, 15, ncols);

                    simdfunc::contract_primitives(buffer, 624, 102, 15, ncols);

                    simdfunc::contract_primitives(buffer, 654, 319, 21, ncols);

                    simdfunc::contract_primitives(buffer, 696, 340, 21, ncols);

                    simdfunc::contract_primitives(buffer, 738, 361, 21, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 385, 382, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 391, 388, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 397, 394, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 403, 400, 3, 1, nmax);

        simdtrf::transform_s_inner(buffer, 412, 406, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 424, 418, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 436, 430, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 448, 442, 6, 1, nmax);

        simdtrf::transform_s_inner(buffer, 464, 454, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 484, 474, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 504, 494, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 524, 514, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 549, 534, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 579, 564, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 609, 594, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 639, 624, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 675, 654, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 717, 696, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 759, 738, 21, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 780, 385, 403, 412,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 789, 391, 403, 424,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 798, 397, 403, 436,
                                                       1, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 807, 403, 448, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 816, 412, 448, 464,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 834, 424, 448, 484,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 852, 436, 448, 504,
                                                       1, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 870, 448, 524, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 888, 464, 524, 549,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 918, 484, 524, 579,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 948, 504, 524, 609,
                                                       1, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 978, 524, 639, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 1008, 549, 639, 675,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 1053, 579, 639, 717,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 1098, 609, 639, 759,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 1143, 780, 807, 816,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 1161, 789, 807, 834,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 1179, 798, 807, 852,
                                                       1, nmax);

        simdtrf::compute_hrr_pd_out_of_first(buffer, coordinates, 1197, 807, 870, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 1215, 816, 870, 888,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 1251, 834, 870, 918,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 1287, 852, 870, 948,
                                                       1, nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 1323, 870, 978, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 1359, 888, 978, 1008,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 1419, 918, 978, 1053,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 1479, 948, 978, 1098,
                                                       1, nmax);

        simdtrf::compute_hrr_geom_100x_pf_out_of_first(buffer, coordinates, 1539, 1143, 1197,
                                                       1215, 1, nmax);

        simdtrf::compute_hrr_geom_100y_pf_out_of_first(buffer, coordinates, 1569, 1161, 1197,
                                                       1251, 1, nmax);

        simdtrf::compute_hrr_geom_100z_pf_out_of_first(buffer, coordinates, 1599, 1179, 1197,
                                                       1287, 1, nmax);

        simdtrf::compute_hrr_pf_out_of_first(buffer, coordinates, 1629, 1197, 1323, 1, nmax);

        simdtrf::compute_hrr_geom_100x_df_out_of_first(buffer, coordinates, 1659, 1215, 1323,
                                                       1359, 1, nmax);

        simdtrf::compute_hrr_geom_100y_df_out_of_first(buffer, coordinates, 1719, 1251, 1323,
                                                       1419, 1, nmax);

        simdtrf::compute_hrr_geom_100z_df_out_of_first(buffer, coordinates, 1779, 1287, 1323,
                                                       1479, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pg_out_of_first(buffer, coordinates, 1839, 1539, 1629,
                                                       1659, 1, nmax);

        simdtrf::compute_hrr_geom_100y_pg_out_of_first(buffer, coordinates, 1884, 1569, 1629,
                                                       1719, 1, nmax);

        simdtrf::compute_hrr_geom_100z_pg_out_of_first(buffer, coordinates, 1929, 1599, 1629,
                                                       1779, 1, nmax);

        simdtrf::transform_g_inner(buffer, 1974, 1839, 3, 1, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 1974, 9, nmax);

        simdtrf::transform_g_inner(buffer, 1974, 1884, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 27 * nvalues + n * npairs, nvalues, buffer, 1974, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 1974, 1929, 3, 1, nmax);

        simdtrf::transform_p_outer(values + 54 * nvalues + n * npairs, nvalues, buffer, 1974, 9,
                                   nmax);
    }

    for (size_t m = 0; m < 81; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
