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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecFGS.hpp"

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
compute_rs_geom_100_fgs_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_fgs_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 9528, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 9528, 1658, 1580, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 8,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 16, 3, 8,
                                                             ncols, fj, i * nprim_b + j, fq);

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

                    simdgeo::geom_f_x(buffer, 998, 74, 278, 1, 1, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 1008, 74, 278, 1, 1, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 1018, 74, 278, 1, 1, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 1028, 116, 353, 1, 1, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 1038, 116, 353, 1, 1, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 1048, 116, 353, 1, 1, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 1058, 158, 428, 1, 1, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 1073, 158, 428, 1, 1, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 1088, 158, 428, 1, 1, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 1103, 218, 512, 1, 1, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 1118, 218, 512, 1, 1, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 1133, 218, 512, 1, 1, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 1148, 278, 596, 1, 1, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 1169, 278, 596, 1, 1, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 1190, 278, 596, 1, 1, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 1211, 353, 680, 1, 1, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 1232, 353, 680, 1, 1, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 1253, 353, 680, 1, 1, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 1274, 428, 764, 1, 1, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 1302, 428, 764, 1, 1, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 1330, 428, 764, 1, 1, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 1358, 512, 836, 1, 1, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 1386, 512, 836, 1, 1, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 1414, 512, 836, 1, 1, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 1442, 596, 908, 1, 1, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 1478, 596, 908, 1, 1, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 1514, 596, 908, 1, 1, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 1550, 680, 953, 1, 1, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 1586, 680, 953, 1, 1, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 1622, 680, 953, 1, 1, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 1658, 998, 10, ncols);

                    simdfunc::contract_primitives(buffer, 1678, 1008, 10, ncols);

                    simdfunc::contract_primitives(buffer, 1698, 1018, 10, ncols);

                    simdfunc::contract_primitives(buffer, 1718, 158, 10, ncols);

                    simdfunc::contract_primitives(buffer, 1738, 1028, 10, ncols);

                    simdfunc::contract_primitives(buffer, 1758, 1038, 10, ncols);

                    simdfunc::contract_primitives(buffer, 1778, 1048, 10, ncols);

                    simdfunc::contract_primitives(buffer, 1798, 218, 10, ncols);

                    simdfunc::contract_primitives(buffer, 1818, 1058, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1848, 1073, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1878, 1088, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1908, 278, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1938, 1103, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1968, 1118, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1998, 1133, 15, ncols);

                    simdfunc::contract_primitives(buffer, 2028, 353, 15, ncols);

                    simdfunc::contract_primitives(buffer, 2058, 1148, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2100, 1169, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2142, 1190, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2184, 428, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2226, 1211, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2268, 1232, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2310, 1253, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2352, 512, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2394, 1274, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2450, 1302, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2506, 1330, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2562, 596, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2618, 1358, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2674, 1386, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2730, 1414, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2786, 680, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2842, 1442, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2914, 1478, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2986, 1514, 36, ncols);

                    simdfunc::contract_primitives(buffer, 3058, 1550, 36, ncols);

                    simdfunc::contract_primitives(buffer, 3130, 1586, 36, ncols);

                    simdfunc::contract_primitives(buffer, 3202, 1622, 36, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 1668, 1658, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1688, 1678, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1708, 1698, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1728, 1718, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1748, 1738, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1768, 1758, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1788, 1778, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1808, 1798, 10, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1833, 1818, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1863, 1848, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1893, 1878, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1923, 1908, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1953, 1938, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1983, 1968, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2013, 1998, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2043, 2028, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2079, 2058, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2121, 2100, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2163, 2142, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2205, 2184, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2247, 2226, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2289, 2268, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2331, 2310, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2373, 2352, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2422, 2394, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2478, 2450, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2534, 2506, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2590, 2562, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2646, 2618, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2702, 2674, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2758, 2730, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2814, 2786, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2878, 2842, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2950, 2914, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3022, 2986, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3094, 3058, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3166, 3130, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3238, 3202, 36, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 3274, 1668, 1728,
                                                       1833, 1, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 3304, 1688, 1728,
                                                       1863, 1, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 3334, 1708, 1728,
                                                       1893, 1, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 3364, 1728, 1923, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 3394, 1748, 1808,
                                                       1953, 1, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 3424, 1768, 1808,
                                                       1983, 1, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 3454, 1788, 1808,
                                                       2013, 1, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 3484, 1808, 2043, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 3514, 1833, 1923,
                                                       2079, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 3559, 1863, 1923,
                                                       2121, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 3604, 1893, 1923,
                                                       2163, 1, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 3649, 1923, 2205, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 3694, 1953, 2043,
                                                       2247, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 3739, 1983, 2043,
                                                       2289, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 3784, 2013, 2043,
                                                       2331, 1, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 3829, 2043, 2373, 1, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 3874, 2079, 2205,
                                                       2422, 1, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 3937, 2121, 2205,
                                                       2478, 1, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 4000, 2163, 2205,
                                                       2534, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 4063, 2205, 2590, 1, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 4126, 2247, 2373,
                                                       2646, 1, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 4189, 2289, 2373,
                                                       2702, 1, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 4252, 2331, 2373,
                                                       2758, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 4315, 2373, 2814, 1, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 4378, 2422, 2590,
                                                       2878, 1, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 4462, 2478, 2590,
                                                       2950, 1, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 4546, 2534, 2590,
                                                       3022, 1, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 4630, 2646, 2814,
                                                       3094, 1, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 4714, 2702, 2814,
                                                       3166, 1, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 4798, 2758, 2814,
                                                       3238, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 4882, 3274, 3364,
                                                       3514, 1, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 4942, 3304, 3364,
                                                       3559, 1, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 5002, 3334, 3364,
                                                       3604, 1, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 5062, 3364, 3649, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 5122, 3394, 3484,
                                                       3694, 1, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 5182, 3424, 3484,
                                                       3739, 1, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 5242, 3454, 3484,
                                                       3784, 1, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 5302, 3484, 3829, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 5362, 3514, 3649,
                                                       3874, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 5452, 3559, 3649,
                                                       3937, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 5542, 3604, 3649,
                                                       4000, 1, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 5632, 3649, 4063, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 5722, 3694, 3829,
                                                       4126, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 5812, 3739, 3829,
                                                       4189, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 5902, 3784, 3829,
                                                       4252, 1, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 5992, 3829, 4315, 1, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 6082, 3874, 4063,
                                                       4378, 1, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 6208, 3937, 4063,
                                                       4462, 1, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 6334, 4000, 4063,
                                                       4546, 1, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 6460, 4126, 4315,
                                                       4630, 1, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 6586, 4189, 4315,
                                                       4714, 1, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 6712, 4252, 4315,
                                                       4798, 1, nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 6838, 4882, 5062,
                                                       5362, 1, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 6938, 4942, 5062,
                                                       5452, 1, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 7038, 5002, 5062,
                                                       5542, 1, nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 7138, 5062, 5632, 1, nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 7238, 5122, 5302,
                                                       5722, 1, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 7338, 5182, 5302,
                                                       5812, 1, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 7438, 5242, 5302,
                                                       5902, 1, nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 7538, 5302, 5992, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 7638, 5362, 5632,
                                                       6082, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 7788, 5452, 5632,
                                                       6208, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 7938, 5542, 5632,
                                                       6334, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 8088, 5722, 5992,
                                                       6460, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 8238, 5812, 5992,
                                                       6586, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 8388, 5902, 5992,
                                                       6712, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fg_out_of_first(buffer, coordinates, 8538, 6838, 7138,
                                                       7638, 1, nmax);

        simdtrf::compute_hrr_geom_100y_fg_out_of_first(buffer, coordinates, 8688, 6938, 7138,
                                                       7788, 1, nmax);

        simdtrf::compute_hrr_geom_100z_fg_out_of_first(buffer, coordinates, 8838, 7038, 7138,
                                                       7938, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fg_out_of_first(buffer, coordinates, 8988, 7238, 7538,
                                                       8088, 1, nmax);

        simdtrf::compute_hrr_geom_100y_fg_out_of_first(buffer, coordinates, 9138, 7338, 7538,
                                                       8238, 1, nmax);

        simdtrf::compute_hrr_geom_100z_fg_out_of_first(buffer, coordinates, 9288, 7438, 7538,
                                                       8388, 1, nmax);

        simdtrf::transform_g_inner(buffer, 9438, 8988, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 9438, 9, nmax);

        simdtrf::transform_g_inner(buffer, 9438, 9138, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 9438, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 9438, 9288, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 126 * nvalues + n * npairs, nvalues, buffer, 9438, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 9438, 8538, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 189 * nvalues + n * npairs, nvalues, buffer, 9438, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 9438, 8688, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 252 * nvalues + n * npairs, nvalues, buffer, 9438, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 9438, 8838, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 9438, 9,
                                   nmax);
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
