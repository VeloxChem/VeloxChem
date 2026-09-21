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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecGFS.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
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
#include "SimdTransformS.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_gfs_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_gfs_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 6903, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 6903, 1598, 1420, dimensions);

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

                    simdgeo::geom_g_x(buffer, 998, 158, 428, 1, 1, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 1013, 158, 428, 1, 1, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 1028, 158, 428, 1, 1, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 1043, 218, 512, 1, 1, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 1058, 218, 512, 1, 1, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 1073, 218, 512, 1, 1, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 1088, 278, 596, 1, 1, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 1109, 278, 596, 1, 1, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 1130, 278, 596, 1, 1, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 1151, 353, 680, 1, 1, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 1172, 353, 680, 1, 1, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 1193, 353, 680, 1, 1, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 1214, 428, 764, 1, 1, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 1242, 428, 764, 1, 1, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 1270, 428, 764, 1, 1, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 1298, 512, 836, 1, 1, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 1326, 512, 836, 1, 1, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 1354, 512, 836, 1, 1, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 1382, 596, 908, 1, 1, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 1418, 596, 908, 1, 1, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 1454, 596, 908, 1, 1, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 1490, 680, 953, 1, 1, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 1526, 680, 953, 1, 1, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 1562, 680, 953, 1, 1, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 1598, 998, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1628, 1013, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1658, 1028, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1688, 278, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1718, 1043, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1748, 1058, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1778, 1073, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1808, 353, 15, ncols);

                    simdfunc::contract_primitives(buffer, 1838, 1088, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1880, 1109, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1922, 1130, 21, ncols);

                    simdfunc::contract_primitives(buffer, 1964, 428, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2006, 1151, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2048, 1172, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2090, 1193, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2132, 512, 21, ncols);

                    simdfunc::contract_primitives(buffer, 2174, 1214, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2230, 1242, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2286, 1270, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2342, 596, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2398, 1298, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2454, 1326, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2510, 1354, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2566, 680, 28, ncols);

                    simdfunc::contract_primitives(buffer, 2622, 1382, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2694, 1418, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2766, 1454, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2838, 1490, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2910, 1526, 36, ncols);

                    simdfunc::contract_primitives(buffer, 2982, 1562, 36, ncols);
                }
            }
        }

        simdtrf::transform_s_inner(buffer, 1613, 1598, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1643, 1628, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1673, 1658, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1703, 1688, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1733, 1718, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1763, 1748, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1793, 1778, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1823, 1808, 15, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1859, 1838, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1901, 1880, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1943, 1922, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 1985, 1964, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2027, 2006, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2069, 2048, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2111, 2090, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2153, 2132, 21, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2202, 2174, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2258, 2230, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2314, 2286, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2370, 2342, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2426, 2398, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2482, 2454, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2538, 2510, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2594, 2566, 28, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2658, 2622, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2730, 2694, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2802, 2766, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2874, 2838, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 2946, 2910, 36, 1, nmax);

        simdtrf::transform_s_inner(buffer, 3018, 2982, 36, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 3054, 1613, 1703,
                                                       1859, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 3099, 1643, 1703,
                                                       1901, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 3144, 1673, 1703,
                                                       1943, 1, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 3189, 1703, 1985, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 3234, 1733, 1823,
                                                       2027, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 3279, 1763, 1823,
                                                       2069, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 3324, 1793, 1823,
                                                       2111, 1, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 3369, 1823, 2153, 1, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 3414, 1859, 1985,
                                                       2202, 1, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 3477, 1901, 1985,
                                                       2258, 1, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 3540, 1943, 1985,
                                                       2314, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 3603, 1985, 2370, 1, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 3666, 2027, 2153,
                                                       2426, 1, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 3729, 2069, 2153,
                                                       2482, 1, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 3792, 2111, 2153,
                                                       2538, 1, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 3855, 2153, 2594, 1, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 3918, 2202, 2370,
                                                       2658, 1, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 4002, 2258, 2370,
                                                       2730, 1, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 4086, 2314, 2370,
                                                       2802, 1, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 4170, 2426, 2594,
                                                       2874, 1, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 4254, 2482, 2594,
                                                       2946, 1, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 4338, 2538, 2594,
                                                       3018, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 4422, 3054, 3189,
                                                       3414, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 4512, 3099, 3189,
                                                       3477, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 4602, 3144, 3189,
                                                       3540, 1, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 4692, 3189, 3603, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 4782, 3234, 3369,
                                                       3666, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 4872, 3279, 3369,
                                                       3729, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 4962, 3324, 3369,
                                                       3792, 1, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 5052, 3369, 3855, 1, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 5142, 3414, 3603,
                                                       3918, 1, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 5268, 3477, 3603,
                                                       4002, 1, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 5394, 3540, 3603,
                                                       4086, 1, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 5520, 3666, 3855,
                                                       4170, 1, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 5646, 3729, 3855,
                                                       4254, 1, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 5772, 3792, 3855,
                                                       4338, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 5898, 4422, 4692,
                                                       5142, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 6048, 4512, 4692,
                                                       5268, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 6198, 4602, 4692,
                                                       5394, 1, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 6348, 4782, 5052,
                                                       5520, 1, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 6498, 4872, 5052,
                                                       5646, 1, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 6648, 4962, 5052,
                                                       5772, 1, nmax);

        simdtrf::transform_f_inner(buffer, 6798, 6348, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 6798, 7, nmax);

        simdtrf::transform_f_inner(buffer, 6798, 6498, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 6798, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 6798, 6648, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 126 * nvalues + n * npairs, nvalues, buffer, 6798, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 6798, 5898, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 189 * nvalues + n * npairs, nvalues, buffer, 6798, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 6798, 6048, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 252 * nvalues + n * npairs, nvalues, buffer, 6798, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 6798, 6198, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 6798, 7,
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
