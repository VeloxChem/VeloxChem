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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecPFP.hpp"

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
#include "SimdGeometryP1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferGeom100XDD.hpp"
#include "SimdTransferGeom100XDP.hpp"
#include "SimdTransferGeom100XFP.hpp"
#include "SimdTransferGeom100XPD.hpp"
#include "SimdTransferGeom100XPF.hpp"
#include "SimdTransferGeom100XPP.hpp"
#include "SimdTransferGeom100YDD.hpp"
#include "SimdTransferGeom100YDP.hpp"
#include "SimdTransferGeom100YFP.hpp"
#include "SimdTransferGeom100YPD.hpp"
#include "SimdTransferGeom100YPF.hpp"
#include "SimdTransferGeom100YPP.hpp"
#include "SimdTransferGeom100ZDD.hpp"
#include "SimdTransferGeom100ZDP.hpp"
#include "SimdTransferGeom100ZFP.hpp"
#include "SimdTransferGeom100ZPD.hpp"
#include "SimdTransferGeom100ZPF.hpp"
#include "SimdTransferGeom100ZPP.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransferPP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_pfp_three_center_electron_repulsion(double               *values,
                                                        const size_t          npairs,
                                                        const size_t          natoms,
                                                        const CBasisFunction &a_function,
                                                        const CBasisFunction &b_function,
                                                        const CBasisFunction &c_function,
                                                        const CSimdMatrix    &coordinates,
                                                        const CSimdMatrix    &c_coordinates,
                                                        CSimdMatrix          &buffer,
                                                        const double          omega,
                                                        const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_rs_geom_100_pfp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 5531, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 378 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 5531, 1208, 1407, dimensions);

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

                    simdfunc::compute_t3c_erf_boys_function(buffer, coordinates, 6, 3, {1, 2, 3,
                                                            4, 5, 6}, ncols, fj, i * nprim_b + j,
                                                            fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 13, 3, {1, 2, 3, 4,
                                                        5, 6}, ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 20, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 23, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 7, 8,
                                                                       20, 23, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 8, 9,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 9, 10,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 10, 11,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 14, 15,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 15, 16,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 16, 17,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 17, 18,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 20, 23,
                                                                       50, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 23, 26,
                                                                       56, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 26, 29,
                                                                       62, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 35, 38,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 38, 41,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 41, 44,
                                                                       86, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 50, 56,
                                                                       98, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 173, 0, 3, 56, 62,
                                                                       108, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 74, 80,
                                                                       128, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 203, 0, 3, 80, 86,
                                                                       138, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 98,
                                                                       108, 158, 173, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 239, 0, 3, 128,
                                                                       138, 188, 203, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 260, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 263, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 266, 3, 7, 20,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 275, 3, 14, 35,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 284, 3, 20, 50,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 302, 3, 35, 74,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 320, 3, 50, 98,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 350, 3, 74, 128,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 380, 3, 98, 158,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 425, 3, 128, 188,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 470, 3, 158, 218,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 533, 3, 188, 239,
                                                                       ncols, p, q);

                    simdgeo::geom_p_x(buffer, 596, 260, 284, 1, 3, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 605, 260, 284, 1, 3, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 614, 260, 284, 1, 3, ncols, alpha);

                    simdgeo::geom_p_x(buffer, 623, 263, 302, 1, 3, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 632, 263, 302, 1, 3, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 641, 263, 302, 1, 3, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 650, 266, 320, 1, 3, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 668, 266, 320, 1, 3, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 686, 266, 320, 1, 3, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 704, 275, 350, 1, 3, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 722, 275, 350, 1, 3, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 740, 275, 350, 1, 3, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 758, 284, 380, 1, 3, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 788, 284, 380, 1, 3, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 818, 284, 380, 1, 3, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 848, 302, 425, 1, 3, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 878, 302, 425, 1, 3, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 908, 302, 425, 1, 3, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 938, 320, 470, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 983, 320, 470, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 1028, 320, 470, 1, 3, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 1073, 350, 533, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 1118, 350, 533, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 1163, 350, 533, 1, 3, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 1208, 596, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1226, 605, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1244, 614, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1262, 266, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1280, 623, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1298, 632, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1316, 641, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1334, 275, 9, ncols);

                    simdfunc::contract_primitives(buffer, 1352, 650, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1388, 668, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1424, 686, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1460, 284, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1496, 704, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1532, 722, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1568, 740, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1604, 302, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1640, 758, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1700, 788, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1760, 818, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1820, 320, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1880, 848, 30, ncols);

                    simdfunc::contract_primitives(buffer, 1940, 878, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2000, 908, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2060, 350, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2120, 938, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2210, 983, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2300, 1028, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2390, 1073, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2480, 1118, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2570, 1163, 45, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 1217, 1208, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1235, 1226, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1253, 1244, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1271, 1262, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1289, 1280, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1307, 1298, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1325, 1316, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1343, 1334, 3, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1370, 1352, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1406, 1388, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1442, 1424, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1478, 1460, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1514, 1496, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1550, 1532, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1586, 1568, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1622, 1604, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1670, 1640, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1730, 1700, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1790, 1760, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1850, 1820, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1910, 1880, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1970, 1940, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2030, 2000, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2090, 2060, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2165, 2120, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2255, 2210, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2345, 2300, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2435, 2390, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2525, 2480, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2615, 2570, 15, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 2660, 1217, 1271,
                                                       1370, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 2687, 1235, 1271,
                                                       1406, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 2714, 1253, 1271,
                                                       1442, 3, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 2741, 1271, 1478, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 2768, 1289, 1343,
                                                       1514, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 2795, 1307, 1343,
                                                       1550, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 2822, 1325, 1343,
                                                       1586, 3, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 2849, 1343, 1622, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 2876, 1370, 1478,
                                                       1670, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 2930, 1406, 1478,
                                                       1730, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 2984, 1442, 1478,
                                                       1790, 3, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 3038, 1478, 1850, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 3092, 1514, 1622,
                                                       1910, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 3146, 1550, 1622,
                                                       1970, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 3200, 1586, 1622,
                                                       2030, 3, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 3254, 1622, 2090, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 3308, 1670, 1850,
                                                       2165, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 3398, 1730, 1850,
                                                       2255, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 3488, 1790, 1850,
                                                       2345, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 3578, 1910, 2090,
                                                       2435, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 3668, 1970, 2090,
                                                       2525, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 3758, 2030, 2090,
                                                       2615, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 3848, 2660, 2741,
                                                       2876, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 3902, 2687, 2741,
                                                       2930, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 3956, 2714, 2741,
                                                       2984, 3, nmax);

        simdtrf::compute_hrr_pd_out_of_first(buffer, coordinates, 4010, 2741, 3038, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 4064, 2768, 2849,
                                                       3092, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 4118, 2795, 2849,
                                                       3146, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 4172, 2822, 2849,
                                                       3200, 3, nmax);

        simdtrf::compute_hrr_pd_out_of_first(buffer, coordinates, 4226, 2849, 3254, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 4280, 2876, 3038,
                                                       3308, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 4388, 2930, 3038,
                                                       3398, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 4496, 2984, 3038,
                                                       3488, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 4604, 3092, 3254,
                                                       3578, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 4712, 3146, 3254,
                                                       3668, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 4820, 3200, 3254,
                                                       3758, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pf_out_of_first(buffer, coordinates, 4928, 3848, 4010,
                                                       4280, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pf_out_of_first(buffer, coordinates, 5018, 3902, 4010,
                                                       4388, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pf_out_of_first(buffer, coordinates, 5108, 3956, 4010,
                                                       4496, 3, nmax);

        simdtrf::compute_hrr_geom_100x_pf_out_of_first(buffer, coordinates, 5198, 4064, 4226,
                                                       4604, 3, nmax);

        simdtrf::compute_hrr_geom_100y_pf_out_of_first(buffer, coordinates, 5288, 4118, 4226,
                                                       4712, 3, nmax);

        simdtrf::compute_hrr_geom_100z_pf_out_of_first(buffer, coordinates, 5378, 4172, 4226,
                                                       4820, 3, nmax);

        simdtrf::transform_f_inner(buffer, 5468, 5198, 3, 3, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 5468, 21, nmax);

        simdtrf::transform_f_inner(buffer, 5468, 5288, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 5468, 21,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 5468, 5378, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 126 * nvalues + n * npairs, nvalues, buffer, 5468,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 5468, 4928, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 189 * nvalues + n * npairs, nvalues, buffer, 5468,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 5468, 5018, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 252 * nvalues + n * npairs, nvalues, buffer, 5468,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 5468, 5108, 3, 3, nmax);

        simdtrf::transform_p_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 5468,
                                   21, nmax);
    }

    for (size_t m = 0; m < 378; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
