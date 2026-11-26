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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecFGP.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
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
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_fgp_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_fgp_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 27554, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 27554, 3944, 4740, dimensions);

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

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 998, 3, 26, 74,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1016, 3, 50, 116,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1034, 3, 74, 158,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1064, 3, 116, 218,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1094, 3, 158, 278,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1139, 3, 218, 353,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1184, 3, 278, 428,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1247, 3, 353, 512,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 1310, 3, 428, 596,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 1394, 3, 512, 680,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 1478, 3, 596, 764,
                                                                       ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 1586, 3, 680, 836,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 1694, 3, 764, 908,
                                                                       ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 1829, 3, 836, 953,
                                                                       ncols, p, q);

                    simdgeo::geom_f_x(buffer, 1964, 998, 1094, 1, 3, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 1994, 998, 1094, 1, 3, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 2024, 998, 1094, 1, 3, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 2054, 1016, 1139, 1, 3, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 2084, 1016, 1139, 1, 3, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 2114, 1016, 1139, 1, 3, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 2144, 1034, 1184, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 2189, 1034, 1184, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 2234, 1034, 1184, 1, 3, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 2279, 1064, 1247, 1, 3, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 2324, 1064, 1247, 1, 3, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 2369, 1064, 1247, 1, 3, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 2414, 1094, 1310, 1, 3, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 2477, 1094, 1310, 1, 3, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 2540, 1094, 1310, 1, 3, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 2603, 1139, 1394, 1, 3, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 2666, 1139, 1394, 1, 3, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 2729, 1139, 1394, 1, 3, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 2792, 1184, 1478, 1, 3, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 2876, 1184, 1478, 1, 3, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 2960, 1184, 1478, 1, 3, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 3044, 1247, 1586, 1, 3, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 3128, 1247, 1586, 1, 3, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 3212, 1247, 1586, 1, 3, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 3296, 1310, 1694, 1, 3, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 3404, 1310, 1694, 1, 3, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 3512, 1310, 1694, 1, 3, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 3620, 1394, 1829, 1, 3, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 3728, 1394, 1829, 1, 3, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 3836, 1394, 1829, 1, 3, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 3944, 1964, 30, ncols);

                    simdfunc::contract_primitives(buffer, 4004, 1994, 30, ncols);

                    simdfunc::contract_primitives(buffer, 4064, 2024, 30, ncols);

                    simdfunc::contract_primitives(buffer, 4124, 1034, 30, ncols);

                    simdfunc::contract_primitives(buffer, 4184, 2054, 30, ncols);

                    simdfunc::contract_primitives(buffer, 4244, 2084, 30, ncols);

                    simdfunc::contract_primitives(buffer, 4304, 2114, 30, ncols);

                    simdfunc::contract_primitives(buffer, 4364, 1064, 30, ncols);

                    simdfunc::contract_primitives(buffer, 4424, 2144, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4514, 2189, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4604, 2234, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4694, 1094, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4784, 2279, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4874, 2324, 45, ncols);

                    simdfunc::contract_primitives(buffer, 4964, 2369, 45, ncols);

                    simdfunc::contract_primitives(buffer, 5054, 1139, 45, ncols);

                    simdfunc::contract_primitives(buffer, 5144, 2414, 63, ncols);

                    simdfunc::contract_primitives(buffer, 5270, 2477, 63, ncols);

                    simdfunc::contract_primitives(buffer, 5396, 2540, 63, ncols);

                    simdfunc::contract_primitives(buffer, 5522, 1184, 63, ncols);

                    simdfunc::contract_primitives(buffer, 5648, 2603, 63, ncols);

                    simdfunc::contract_primitives(buffer, 5774, 2666, 63, ncols);

                    simdfunc::contract_primitives(buffer, 5900, 2729, 63, ncols);

                    simdfunc::contract_primitives(buffer, 6026, 1247, 63, ncols);

                    simdfunc::contract_primitives(buffer, 6152, 2792, 84, ncols);

                    simdfunc::contract_primitives(buffer, 6320, 2876, 84, ncols);

                    simdfunc::contract_primitives(buffer, 6488, 2960, 84, ncols);

                    simdfunc::contract_primitives(buffer, 6656, 1310, 84, ncols);

                    simdfunc::contract_primitives(buffer, 6824, 3044, 84, ncols);

                    simdfunc::contract_primitives(buffer, 6992, 3128, 84, ncols);

                    simdfunc::contract_primitives(buffer, 7160, 3212, 84, ncols);

                    simdfunc::contract_primitives(buffer, 7328, 1394, 84, ncols);

                    simdfunc::contract_primitives(buffer, 7496, 3296, 108, ncols);

                    simdfunc::contract_primitives(buffer, 7712, 3404, 108, ncols);

                    simdfunc::contract_primitives(buffer, 7928, 3512, 108, ncols);

                    simdfunc::contract_primitives(buffer, 8144, 3620, 108, ncols);

                    simdfunc::contract_primitives(buffer, 8360, 3728, 108, ncols);

                    simdfunc::contract_primitives(buffer, 8576, 3836, 108, ncols);
                }
            }
        }

        simdtrf::transform_p_inner(buffer, 3974, 3944, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4034, 4004, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4094, 4064, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4154, 4124, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4214, 4184, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4274, 4244, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4334, 4304, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4394, 4364, 10, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4469, 4424, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4559, 4514, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4649, 4604, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4739, 4694, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4829, 4784, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 4919, 4874, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5009, 4964, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5099, 5054, 15, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5207, 5144, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5333, 5270, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5459, 5396, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5585, 5522, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5711, 5648, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5837, 5774, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 5963, 5900, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6089, 6026, 21, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6236, 6152, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6404, 6320, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6572, 6488, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6740, 6656, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 6908, 6824, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7076, 6992, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7244, 7160, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7412, 7328, 28, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7604, 7496, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 7820, 7712, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 8036, 7928, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 8252, 8144, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 8468, 8360, 36, 1, nmax);

        simdtrf::transform_p_inner(buffer, 8684, 8576, 36, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 8792, 3974, 4154,
                                                       4469, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 8882, 4034, 4154,
                                                       4559, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 8972, 4094, 4154,
                                                       4649, 3, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 9062, 4154, 4739, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 9152, 4214, 4394,
                                                       4829, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 9242, 4274, 4394,
                                                       4919, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 9332, 4334, 4394,
                                                       5009, 3, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 9422, 4394, 5099, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 9512, 4469, 4739,
                                                       5207, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 9647, 4559, 4739,
                                                       5333, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 9782, 4649, 4739,
                                                       5459, 3, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 9917, 4739, 5585, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 10052, 4829, 5099,
                                                       5711, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 10187, 4919, 5099,
                                                       5837, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 10322, 5009, 5099,
                                                       5963, 3, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 10457, 5099, 6089, 3, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 10592, 5207, 5585,
                                                       6236, 3, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 10781, 5333, 5585,
                                                       6404, 3, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 10970, 5459, 5585,
                                                       6572, 3, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 11159, 5585, 6740, 3, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 11348, 5711, 6089,
                                                       6908, 3, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 11537, 5837, 6089,
                                                       7076, 3, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 11726, 5963, 6089,
                                                       7244, 3, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 11915, 6089, 7412, 3, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 12104, 6236, 6740,
                                                       7604, 3, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 12356, 6404, 6740,
                                                       7820, 3, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 12608, 6572, 6740,
                                                       8036, 3, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 12860, 6908, 7412,
                                                       8252, 3, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 13112, 7076, 7412,
                                                       8468, 3, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 13364, 7244, 7412,
                                                       8684, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 13616, 8792, 9062,
                                                       9512, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 13796, 8882, 9062,
                                                       9647, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 13976, 8972, 9062,
                                                       9782, 3, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 14156, 9062, 9917, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 14336, 9152, 9422,
                                                       10052, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 14516, 9242, 9422,
                                                       10187, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 14696, 9332, 9422,
                                                       10322, 3, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 14876, 9422, 10457, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 15056, 9512, 9917,
                                                       10592, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 15326, 9647, 9917,
                                                       10781, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 15596, 9782, 9917,
                                                       10970, 3, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 15866, 9917, 11159, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 16136, 10052, 10457,
                                                       11348, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 16406, 10187, 10457,
                                                       11537, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 16676, 10322, 10457,
                                                       11726, 3, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 16946, 10457, 11915, 3, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 17216, 10592, 11159,
                                                       12104, 3, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 17594, 10781, 11159,
                                                       12356, 3, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 17972, 10970, 11159,
                                                       12608, 3, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 18350, 11348, 11915,
                                                       12860, 3, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 18728, 11537, 11915,
                                                       13112, 3, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 19106, 11726, 11915,
                                                       13364, 3, nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 19484, 13616, 14156,
                                                       15056, 3, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 19784, 13796, 14156,
                                                       15326, 3, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 20084, 13976, 14156,
                                                       15596, 3, nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 20384, 14156, 15866, 3, nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 20684, 14336, 14876,
                                                       16136, 3, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 20984, 14516, 14876,
                                                       16406, 3, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 21284, 14696, 14876,
                                                       16676, 3, nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 21584, 14876, 16946, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 21884, 15056, 15866,
                                                       17216, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 22334, 15326, 15866,
                                                       17594, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 22784, 15596, 15866,
                                                       17972, 3, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 23234, 16136, 16946,
                                                       18350, 3, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 23684, 16406, 16946,
                                                       18728, 3, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 24134, 16676, 16946,
                                                       19106, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fg_out_of_first(buffer, coordinates, 24584, 19484, 20384,
                                                       21884, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fg_out_of_first(buffer, coordinates, 25034, 19784, 20384,
                                                       22334, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fg_out_of_first(buffer, coordinates, 25484, 20084, 20384,
                                                       22784, 3, nmax);

        simdtrf::compute_hrr_geom_100x_fg_out_of_first(buffer, coordinates, 25934, 20684, 21584,
                                                       23234, 3, nmax);

        simdtrf::compute_hrr_geom_100y_fg_out_of_first(buffer, coordinates, 26384, 20984, 21584,
                                                       23684, 3, nmax);

        simdtrf::compute_hrr_geom_100z_fg_out_of_first(buffer, coordinates, 26834, 21284, 21584,
                                                       24134, 3, nmax);

        simdtrf::transform_g_inner(buffer, 27284, 25934, 10, 3, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 27284, 27, nmax);

        simdtrf::transform_g_inner(buffer, 27284, 26384, 10, 3, nmax);

        simdtrf::transform_f_outer(values + 189 * nvalues + n * npairs, nvalues, buffer, 27284,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 27284, 26834, 10, 3, nmax);

        simdtrf::transform_f_outer(values + 378 * nvalues + n * npairs, nvalues, buffer, 27284,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 27284, 24584, 10, 3, nmax);

        simdtrf::transform_f_outer(values + 567 * nvalues + n * npairs, nvalues, buffer, 27284,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 27284, 25034, 10, 3, nmax);

        simdtrf::transform_f_outer(values + 756 * nvalues + n * npairs, nvalues, buffer, 27284,
                                   27, nmax);

        simdtrf::transform_g_inner(buffer, 27284, 25484, 10, 3, nmax);

        simdtrf::transform_f_outer(values + 945 * nvalues + n * npairs, nvalues, buffer, 27284,
                                   27, nmax);
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
