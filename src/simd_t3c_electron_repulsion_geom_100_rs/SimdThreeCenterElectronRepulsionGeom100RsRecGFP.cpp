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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecGFP.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdTransferGD.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferGeom100XGD.hpp"
#include "SimdTransferGeom100XGF.hpp"
#include "SimdTransferGeom100XGP.hpp"
#include "SimdTransferGeom100XHD.hpp"
#include "SimdTransferGeom100XHP.hpp"
#include "SimdTransferGeom100XIP.hpp"
#include "SimdTransferGeom100YGD.hpp"
#include "SimdTransferGeom100YGF.hpp"
#include "SimdTransferGeom100YGP.hpp"
#include "SimdTransferGeom100YHD.hpp"
#include "SimdTransferGeom100YHP.hpp"
#include "SimdTransferGeom100YIP.hpp"
#include "SimdTransferGeom100ZGD.hpp"
#include "SimdTransferGeom100ZGF.hpp"
#include "SimdTransferGeom100ZGP.hpp"
#include "SimdTransferGeom100ZHD.hpp"
#include "SimdTransferGeom100ZHP.hpp"
#include "SimdTransferGeom100ZIP.hpp"
#include "SimdTransferHP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_gfp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_gfp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 19643, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1134 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 19643, 3728, 4260, dimensions);

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
                                                            4, 5, 6, 7, 8, 9}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 16, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 26, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 29, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 32, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 35, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 38, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 41, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 44, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 47, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 50, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 56, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 59, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 62, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 65, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 68, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 71, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 74, 0, 3, 7, 8,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 80, 0, 3, 8, 9,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 86, 0, 3, 9, 10,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 92, 0, 3, 10, 11,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 98, 0, 3, 11, 12,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 104, 0, 3, 12, 13,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 110, 0, 3, 13, 14,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 116, 0, 3, 17, 18,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 122, 0, 3, 18, 19,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 128, 0, 3, 19, 20,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 134, 0, 3, 20, 21,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 140, 0, 3, 21, 22,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 146, 0, 3, 22, 23,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 152, 0, 3, 23, 24,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 26, 29,
                                                                       74, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 168, 0, 3, 29, 32,
                                                                       80, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 178, 0, 3, 32, 35,
                                                                       86, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 35, 38,
                                                                       92, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 198, 0, 3, 38, 41,
                                                                       98, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 208, 0, 3, 41, 44,
                                                                       104, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 50, 53,
                                                                       116, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 228, 0, 3, 53, 56,
                                                                       122, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 56, 59,
                                                                       128, 134, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 59, 62,
                                                                       134, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 62, 65,
                                                                       140, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 65, 68,
                                                                       146, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 74, 80,
                                                                       158, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 293, 0, 3, 80, 86,
                                                                       168, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 86, 92,
                                                                       178, 188, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 323, 0, 3, 92, 98,
                                                                       188, 198, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 98,
                                                                       104, 198, 208, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 353, 0, 3, 116,
                                                                       122, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 122,
                                                                       128, 228, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 383, 0, 3, 128,
                                                                       134, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 134,
                                                                       140, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 413, 0, 3, 140,
                                                                       146, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 158,
                                                                       168, 278, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 449, 0, 3, 168,
                                                                       178, 293, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 470, 0, 3, 178,
                                                                       188, 308, 323, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 491, 0, 3, 188,
                                                                       198, 323, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 512, 0, 3, 218,
                                                                       228, 353, 368, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 533, 0, 3, 228,
                                                                       238, 368, 383, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 554, 0, 3, 238,
                                                                       248, 383, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 575, 0, 3, 248,
                                                                       258, 398, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 596, 0, 3, 278,
                                                                       293, 428, 449, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 624, 0, 3, 293,
                                                                       308, 449, 470, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 652, 0, 3, 308,
                                                                       323, 470, 491, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 680, 0, 3, 353,
                                                                       368, 512, 533, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 708, 0, 3, 368,
                                                                       383, 533, 554, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 736, 0, 3, 383,
                                                                       398, 554, 575, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 764, 0, 3, 428,
                                                                       449, 596, 624, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 800, 0, 3, 449,
                                                                       470, 624, 652, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 836, 0, 3, 512,
                                                                       533, 680, 708, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 872, 0, 3, 533,
                                                                       554, 708, 736, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 908, 0, 3, 596,
                                                                       624, 764, 800, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 953, 0, 3, 680,
                                                                       708, 836, 872, ncols,
                                                                       gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 998, 3, 74, 158,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1028, 3, 116, 218,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1058, 3, 158, 278,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1103, 3, 218, 353,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1148, 3, 278, 428,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1211, 3, 353, 512,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 1274, 3, 428, 596,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 1358, 3, 512, 680,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 1442, 3, 596, 764,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 1550, 3, 680, 836,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 1658, 3, 764, 908,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 1793, 3, 836, 953,
                                                                       ncols, p, q);

                    simdgeo::geom_g_x(buffer, 1928, 998, 1148, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 1973, 998, 1148, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 2018, 998, 1148, 1, 3, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 2063, 1028, 1211, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 2108, 1028, 1211, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 2153, 1028, 1211, 1, 3, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 2198, 1058, 1274, 1, 3, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 2261, 1058, 1274, 1, 3, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 2324, 1058, 1274, 1, 3, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 2387, 1103, 1358, 1, 3, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 2450, 1103, 1358, 1, 3, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 2513, 1103, 1358, 1, 3, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 2576, 1148, 1442, 1, 3, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 2660, 1148, 1442, 1, 3, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 2744, 1148, 1442, 1, 3, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 2828, 1211, 1550, 1, 3, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 2912, 1211, 1550, 1, 3, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 2996, 1211, 1550, 1, 3, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 3080, 1274, 1658, 1, 3, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 3188, 1274, 1658, 1, 3, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 3296, 1274, 1658, 1, 3, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 3404, 1358, 1793, 1, 3, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 3512, 1358, 1793, 1, 3, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 3620, 1358, 1793, 1, 3, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 3728, 1928, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3818, 1973, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3908, 2018, 45, ncols);

                    simdfunc::contract_primitives(buffer, 3998, 1058, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4088, 2063, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4178, 2108, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4268, 2153, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4358, 1103, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4448, 2198, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4574, 2261, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4700, 2324, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4826, 1148, 63, ncols);

                    simdfunc::contract_primitives(buffer, 4952, 2387, 63, ncols);

                    simdfunc::contract_primitives(buffer, 5078, 2450, 63, ncols);

                    simdfunc::contract_primitives(buffer, 5204, 2513, 63, ncols);

                    simdfunc::contract_primitives(buffer, 5330, 1211, 63, ncols);

                    simdfunc::contract_primitives(buffer, 5456, 2576, 84, ncols);

                    simdfunc::contract_primitives(buffer, 5624, 2660, 84, ncols);

                    simdfunc::contract_primitives(buffer, 5792, 2744, 84, ncols);

                    simdfunc::contract_primitives(buffer, 5960, 1274, 84, ncols);

                    simdfunc::contract_primitives(buffer, 6128, 2828, 84, ncols);

                    simdfunc::contract_primitives(buffer, 6296, 2912, 84, ncols);

                    simdfunc::contract_primitives(buffer, 6464, 2996, 84, ncols);

                    simdfunc::contract_primitives(buffer, 6632, 1358, 84, ncols);

                    simdfunc::contract_primitives(buffer, 6800, 3080, 108, ncols);

                    simdfunc::contract_primitives(buffer, 7016, 3188, 108, ncols);

                    simdfunc::contract_primitives(buffer, 7232, 3296, 108, ncols);

                    simdfunc::contract_primitives(buffer, 7448, 3404, 108, ncols);

                    simdfunc::contract_primitives(buffer, 7664, 3512, 108, ncols);

                    simdfunc::contract_primitives(buffer, 7880, 3620, 108, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 3773, 3728, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3863, 3818, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 3953, 3908, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4043, 3998, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4133, 4088, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4223, 4178, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4313, 4268, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4403, 4358, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4511, 4448, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4637, 4574, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4763, 4700, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4889, 4826, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5015, 4952, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5141, 5078, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5267, 5204, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5393, 5330, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5540, 5456, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5708, 5624, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5876, 5792, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6044, 5960, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6212, 6128, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6380, 6296, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6548, 6464, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6716, 6632, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6908, 6800, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7124, 7016, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7340, 7232, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7556, 7448, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7772, 7664, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7988, 7880, 36, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 8096, 3773, 4043,
                                                       4511, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 8231, 3863, 4043,
                                                       4637, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 8366, 3953, 4043,
                                                       4763, 3, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 8501, 4043, 4889, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 8636, 4133, 4403,
                                                       5015, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 8771, 4223, 4403,
                                                       5141, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 8906, 4313, 4403,
                                                       5267, 3, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 9041, 4403, 5393, 3, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 9176, 4511, 4889,
                                                       5540, 3, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 9365, 4637, 4889,
                                                       5708, 3, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 9554, 4763, 4889,
                                                       5876, 3, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 9743, 4889, 6044, 3, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 9932, 5015, 5393,
                                                       6212, 3, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 10121, 5141, 5393,
                                                       6380, 3, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 10310, 5267, 5393,
                                                       6548, 3, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 10499, 5393, 6716, 3, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 10688, 5540, 6044,
                                                       6908, 3, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 10940, 5708, 6044,
                                                       7124, 3, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 11192, 5876, 6044,
                                                       7340, 3, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 11444, 6212, 6716,
                                                       7556, 3, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 11696, 6380, 6716,
                                                       7772, 3, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 11948, 6548, 6716,
                                                       7988, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 12200, 8096, 8501,
                                                       9176, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 12470, 8231, 8501,
                                                       9365, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 12740, 8366, 8501,
                                                       9554, 3, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 13010, 8501, 9743, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 13280, 8636, 9041,
                                                       9932, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 13550, 8771, 9041,
                                                       10121, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 13820, 8906, 9041,
                                                       10310, 3, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 14090, 9041, 10499, 3, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 14360, 9176, 9743,
                                                       10688, 3, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 14738, 9365, 9743,
                                                       10940, 3, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 15116, 9554, 9743,
                                                       11192, 3, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 15494, 9932, 10499,
                                                       11444, 3, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 15872, 10121, 10499,
                                                       11696, 3, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 16250, 10310, 10499,
                                                       11948, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 16628, 12200, 13010,
                                                       14360, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 17078, 12470, 13010,
                                                       14738, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 17528, 12740, 13010,
                                                       15116, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 17978, 13280, 14090,
                                                       15494, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 18428, 13550, 14090,
                                                       15872, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 18878, 13820, 14090,
                                                       16250, 3, nmax);

        simdtrf::transform_f_inner(buffer, 19328, 17978, 15, 3, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 19328, 21, nmax);

        simdtrf::transform_f_inner(buffer, 19328, 18428, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 189 * nvalues + n * npairs, nvalues, buffer, 19328,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 19328, 18878, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 378 * nvalues + n * npairs, nvalues, buffer, 19328,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 19328, 16628, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 567 * nvalues + n * npairs, nvalues, buffer, 19328,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 19328, 17078, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 756 * nvalues + n * npairs, nvalues, buffer, 19328,
                                   21, nmax);

        simdtrf::transform_f_inner(buffer, 19328, 17528, 15, 3, nmax);

        simdtrf::transform_g_outer(values + 945 * nvalues + n * npairs, nvalues, buffer, 19328,
                                   21, nmax);
    }

    for (size_t m = 0; m < 1134; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
