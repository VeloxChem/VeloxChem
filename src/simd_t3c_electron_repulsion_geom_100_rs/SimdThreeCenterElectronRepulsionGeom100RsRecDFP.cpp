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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecDFP.hpp"

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
#include "SimdTransferDD.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGeom100XDD.hpp"
#include "SimdTransferGeom100XDF.hpp"
#include "SimdTransferGeom100XDP.hpp"
#include "SimdTransferGeom100XFD.hpp"
#include "SimdTransferGeom100XFP.hpp"
#include "SimdTransferGeom100XGP.hpp"
#include "SimdTransferGeom100YDD.hpp"
#include "SimdTransferGeom100YDF.hpp"
#include "SimdTransferGeom100YDP.hpp"
#include "SimdTransferGeom100YFD.hpp"
#include "SimdTransferGeom100YFP.hpp"
#include "SimdTransferGeom100YGP.hpp"
#include "SimdTransferGeom100ZDD.hpp"
#include "SimdTransferGeom100ZDF.hpp"
#include "SimdTransferGeom100ZDP.hpp"
#include "SimdTransferGeom100ZFD.hpp"
#include "SimdTransferGeom100ZFP.hpp"
#include "SimdTransferGeom100ZGP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_dfp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_dfp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 9218, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 630 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 9218, 1862, 2181, dimensions);

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
                                                            4, 5, 6, 7}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 14, 3, {1, 2, 3, 4,
                                                        5, 6, 7}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 22, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 25, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 28, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 31, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 34, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 37, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 40, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 43, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 46, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 49, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 52, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 58, 0, 3, 7, 8,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 64, 0, 3, 8, 9,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 70, 0, 3, 9, 10,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 76, 0, 3, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 82, 0, 3, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 88, 0, 3, 15, 16,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 94, 0, 3, 16, 17,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 100, 0, 3, 17, 18,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 106, 0, 3, 18, 19,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 112, 0, 3, 19, 20,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 118, 0, 3, 22, 25,
                                                                       58, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 25, 28,
                                                                       64, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 28, 31,
                                                                       70, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 31, 34,
                                                                       76, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 40, 43,
                                                                       88, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 168, 0, 3, 43, 46,
                                                                       94, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 178, 0, 3, 46, 49,
                                                                       100, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 49, 52,
                                                                       106, 112, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 198, 0, 3, 58, 64,
                                                                       118, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 213, 0, 3, 64, 70,
                                                                       128, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 228, 0, 3, 70, 76,
                                                                       138, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 243, 0, 3, 88, 94,
                                                                       158, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 94,
                                                                       100, 168, 178, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 273, 0, 3, 100,
                                                                       106, 178, 188, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 118,
                                                                       128, 198, 213, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 309, 0, 3, 128,
                                                                       138, 213, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 330, 0, 3, 158,
                                                                       168, 243, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 351, 0, 3, 168,
                                                                       178, 258, 273, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 372, 0, 3, 198,
                                                                       213, 288, 309, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 400, 0, 3, 243,
                                                                       258, 330, 351, ncols,
                                                                       gamma, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 428, 3, 7, 22,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 437, 3, 15, 40,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 446, 3, 22, 58,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 464, 3, 40, 88,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 482, 3, 58, 118,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 512, 3, 88, 158,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 542, 3, 118, 198,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 587, 3, 158, 243,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 632, 3, 198, 288,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 695, 3, 243, 330,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 758, 3, 288, 372,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 842, 3, 330, 400,
                                                                       ncols, p, q);

                    simdgeo::geom_d_x(buffer, 926, 428, 482, 1, 3, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 944, 428, 482, 1, 3, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 962, 428, 482, 1, 3, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 980, 437, 512, 1, 3, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 998, 437, 512, 1, 3, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 1016, 437, 512, 1, 3, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 1034, 446, 542, 1, 3, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 1064, 446, 542, 1, 3, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 1094, 446, 542, 1, 3, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 1124, 464, 587, 1, 3, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 1154, 464, 587, 1, 3, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 1184, 464, 587, 1, 3, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 1214, 482, 632, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 1259, 482, 632, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 1304, 482, 632, 1, 3, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 1349, 512, 695, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 1394, 512, 695, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 1439, 512, 695, 1, 3, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 1484, 542, 758, 1, 3, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 1547, 542, 758, 1, 3, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 1610, 542, 758, 1, 3, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 1673, 587, 842, 1, 3, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 1736, 587, 842, 1, 3, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 1799, 587, 842, 1, 3, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 1862, 926, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1898, 944, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1934, 962, 18, ncols);

                    simdfunc::contract_primitives(buffer, 1970, 446, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2006, 980, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2042, 998, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2078, 1016, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2114, 464, 18, ncols);

                    simdfunc::contract_primitives(buffer, 2150, 1034, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2210, 1064, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2270, 1094, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2330, 482, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2390, 1124, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2450, 1154, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2510, 1184, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2570, 512, 30, ncols);

                    simdfunc::contract_primitives(buffer, 2630, 1214, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2720, 1259, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2810, 1304, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2900, 542, 45, ncols);

                    simdfunc::contract_primitives(buffer, 2990, 1349, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3080, 1394, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3170, 1439, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3260, 587, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3350, 1484, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3476, 1547, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3602, 1610, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3728, 1673, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3854, 1736, 63, ncols);

                    simdfunc::contract_primitives(buffer, 3980, 1799, 63, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 1880, 1862, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1916, 1898, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1952, 1934, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 1988, 1970, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2024, 2006, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2060, 2042, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2096, 2078, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2132, 2114, 6, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2180, 2150, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2240, 2210, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2300, 2270, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2360, 2330, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2420, 2390, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2480, 2450, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2540, 2510, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2600, 2570, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2675, 2630, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2765, 2720, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2855, 2810, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 2945, 2900, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3035, 2990, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3125, 3080, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3215, 3170, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3305, 3260, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3413, 3350, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3539, 3476, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3665, 3602, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3791, 3728, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3917, 3854, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4043, 3980, 21, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 4106, 1880, 1988,
                                                       2180, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 4160, 1916, 1988,
                                                       2240, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 4214, 1952, 1988,
                                                       2300, 3, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 4268, 1988, 2360, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 4322, 2024, 2132,
                                                       2420, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 4376, 2060, 2132,
                                                       2480, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 4430, 2096, 2132,
                                                       2540, 3, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 4484, 2132, 2600, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 4538, 2180, 2360,
                                                       2675, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 4628, 2240, 2360,
                                                       2765, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 4718, 2300, 2360,
                                                       2855, 3, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 4808, 2360, 2945, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 4898, 2420, 2600,
                                                       3035, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 4988, 2480, 2600,
                                                       3125, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 5078, 2540, 2600,
                                                       3215, 3, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 5168, 2600, 3305, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 5258, 2675, 2945,
                                                       3413, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 5393, 2765, 2945,
                                                       3539, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 5528, 2855, 2945,
                                                       3665, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 5663, 3035, 3305,
                                                       3791, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 5798, 3125, 3305,
                                                       3917, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 5933, 3215, 3305,
                                                       4043, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 6068, 4106, 4268,
                                                       4538, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 6176, 4160, 4268,
                                                       4628, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 6284, 4214, 4268,
                                                       4718, 3, nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 6392, 4268, 4808, 3, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 6500, 4322, 4484,
                                                       4898, 3, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 6608, 4376, 4484,
                                                       4988, 3, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 6716, 4430, 4484,
                                                       5078, 3, nmax);

        simdtrf::compute_hrr_dd_out_of_first(buffer, coordinates, 6824, 4484, 5168, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 6932, 4538, 4808,
                                                       5258, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 7112, 4628, 4808,
                                                       5393, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 7292, 4718, 4808,
                                                       5528, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 7472, 4898, 5168,
                                                       5663, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 7652, 4988, 5168,
                                                       5798, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 7832, 5078, 5168,
                                                       5933, 3, nmax);

        simdtrf::compute_hrr_geom_100x_df_out_of_first(buffer, coordinates, 8012, 6068, 6392,
                                                       6932, 3, nmax);

        simdtrf::compute_hrr_geom_100y_df_out_of_first(buffer, coordinates, 8192, 6176, 6392,
                                                       7112, 3, nmax);

        simdtrf::compute_hrr_geom_100z_df_out_of_first(buffer, coordinates, 8372, 6284, 6392,
                                                       7292, 3, nmax);

        simdtrf::compute_hrr_geom_100x_df_out_of_first(buffer, coordinates, 8552, 6500, 6824,
                                                       7472, 3, nmax);

        simdtrf::compute_hrr_geom_100y_df_out_of_first(buffer, coordinates, 8732, 6608, 6824,
                                                       7652, 3, nmax);

        simdtrf::compute_hrr_geom_100z_df_out_of_first(buffer, coordinates, 8912, 6716, 6824,
                                                       7832, 3, nmax);

        simdtrf::transform_f_inner(buffer, 9092, 8552, 6, 3, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 9092, 21, nmax);

        simdtrf::transform_f_inner(buffer, 9092, 8732, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 105 * nvalues + n * npairs, nvalues, buffer, 9092,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 9092, 8912, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 210 * nvalues + n * npairs, nvalues, buffer, 9092,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 9092, 8012, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 9092,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 9092, 8192, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 420 * nvalues + n * npairs, nvalues, buffer, 9092,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 9092, 8372, 6, 3, nmax);

        simdtrf::transform_d_outer(values + 525 * nvalues + n * npairs, nvalues, buffer, 9092,
                                   21, nmax);
    }

    for (size_t m = 0; m < 630; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
