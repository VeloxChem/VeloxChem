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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecFGD.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecLSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
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
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_fgd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_fgd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 54686, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1890 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 54686, 14528, 8708, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 10,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 18, 3, 10,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 60, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 63, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 66, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 69, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 75, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 78, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 81, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 84, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 87, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 90, 0, 3, 7, 8,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 96, 0, 3, 8, 9,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 9, 10,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 10, 11,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 114, 0, 3, 11, 12,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 120, 0, 3, 12, 13,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 126, 0, 3, 13, 14,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 14, 15,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 15, 16,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 144, 0, 3, 19, 20,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 150, 0, 3, 20, 21,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 156, 0, 3, 21, 22,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 162, 0, 3, 22, 23,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 168, 0, 3, 23, 24,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 174, 0, 3, 24, 25,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 180, 0, 3, 25, 26,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 186, 0, 3, 26, 27,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 27, 28,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 198, 0, 3, 30, 33,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 208, 0, 3, 33, 36,
                                                                       96, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 36, 39,
                                                                       102, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 228, 0, 3, 39, 42,
                                                                       108, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 42, 45,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 248, 0, 3, 45, 48,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 48, 51,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 51, 54,
                                                                       132, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 60, 63,
                                                                       144, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 63, 66,
                                                                       150, 156, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 66, 69,
                                                                       156, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 69, 72,
                                                                       162, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 72, 75,
                                                                       168, 174, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 75, 78,
                                                                       174, 180, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 78, 81,
                                                                       180, 186, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 81, 84,
                                                                       186, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 90, 96,
                                                                       198, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 373, 0, 3, 96,
                                                                       102, 208, 218, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 102,
                                                                       108, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 403, 0, 3, 108,
                                                                       114, 228, 238, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 114,
                                                                       120, 238, 248, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 433, 0, 3, 120,
                                                                       126, 248, 258, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 126,
                                                                       132, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 463, 0, 3, 144,
                                                                       150, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 478, 0, 3, 150,
                                                                       156, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 493, 0, 3, 156,
                                                                       162, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 508, 0, 3, 162,
                                                                       168, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 523, 0, 3, 168,
                                                                       174, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 538, 0, 3, 174,
                                                                       180, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 553, 0, 3, 180,
                                                                       186, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 568, 0, 3, 198,
                                                                       208, 358, 373, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 589, 0, 3, 208,
                                                                       218, 373, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 610, 0, 3, 218,
                                                                       228, 388, 403, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 631, 0, 3, 228,
                                                                       238, 403, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 652, 0, 3, 238,
                                                                       248, 418, 433, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 673, 0, 3, 248,
                                                                       258, 433, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 694, 0, 3, 278,
                                                                       288, 463, 478, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 715, 0, 3, 288,
                                                                       298, 478, 493, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 736, 0, 3, 298,
                                                                       308, 493, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 757, 0, 3, 308,
                                                                       318, 508, 523, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 778, 0, 3, 318,
                                                                       328, 523, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 799, 0, 3, 328,
                                                                       338, 538, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 820, 0, 3, 358,
                                                                       373, 568, 589, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 848, 0, 3, 373,
                                                                       388, 589, 610, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 876, 0, 3, 388,
                                                                       403, 610, 631, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 904, 0, 3, 403,
                                                                       418, 631, 652, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 932, 0, 3, 418,
                                                                       433, 652, 673, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 960, 0, 3, 463,
                                                                       478, 694, 715, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 988, 0, 3, 478,
                                                                       493, 715, 736, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1016, 0, 3, 493,
                                                                       508, 736, 757, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1044, 0, 3, 508,
                                                                       523, 757, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 523,
                                                                       538, 778, 799, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1100, 0, 3, 568,
                                                                       589, 820, 848, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1136, 0, 3, 589,
                                                                       610, 848, 876, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1172, 0, 3, 610,
                                                                       631, 876, 904, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1208, 0, 3, 631,
                                                                       652, 904, 932, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1244, 0, 3, 694,
                                                                       715, 960, 988, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1280, 0, 3, 715,
                                                                       736, 988, 1016, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1316, 0, 3, 736,
                                                                       757, 1016, 1044, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 757,
                                                                       778, 1044, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1388, 0, 3, 820,
                                                                       848, 1100, 1136, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1433, 0, 3, 848,
                                                                       876, 1136, 1172, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1478, 0, 3, 876,
                                                                       904, 1172, 1208, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1523, 0, 3, 960,
                                                                       988, 1244, 1280, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1568, 0, 3, 988,
                                                                       1016, 1280, 1316, ncols,
                                                                       gamma, p, q);

                    compute_prim_lss_three_center_electron_repulsion_0(buffer, 1613, 0, 3, 1016,
                                                                       1044, 1316, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1658, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1661, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1664, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1667, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1670, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1673, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1676, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1679, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1682, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1685, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1688, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1691, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1694, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1697, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1700, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1703, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1706, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1709, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1712, 3, 9, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1721, 3, 10, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1730, 3, 11, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1739, 3, 12, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1748, 3, 13, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1757, 3, 14, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1766, 3, 15, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1775, 3, 16, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1784, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1793, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1802, 3, 23, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1811, 3, 24, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1820, 3, 25, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1829, 3, 26, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1838, 3, 27, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1847, 3, 28, 87,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1856, 3, 36, 102,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1874, 3, 39, 108,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1892, 3, 42, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1910, 3, 45, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1928, 3, 48, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1946, 3, 51, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1964, 3, 54, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1982, 3, 66, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2000, 3, 69, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2018, 3, 72, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2036, 3, 75, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2054, 3, 78, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2072, 3, 81, 186,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2090, 3, 84, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2108, 3, 102, 218,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2138, 3, 108, 228,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2168, 3, 114, 238,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2198, 3, 120, 248,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2228, 3, 126, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2258, 3, 132, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2288, 3, 156, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2318, 3, 162, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2348, 3, 168, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2378, 3, 174, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2408, 3, 180, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2438, 3, 186, 348,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2468, 3, 218, 388,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2513, 3, 228, 403,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2558, 3, 238, 418,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2603, 3, 248, 433,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2648, 3, 258, 448,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2693, 3, 298, 493,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2738, 3, 308, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2783, 3, 318, 523,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2828, 3, 328, 538,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2873, 3, 338, 553,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2918, 3, 388, 610,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2981, 3, 403, 631,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3044, 3, 418, 652,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3107, 3, 433, 673,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3170, 3, 493, 736,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3233, 3, 508, 757,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3296, 3, 523, 778,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3359, 3, 538, 799,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3422, 3, 610, 876,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3506, 3, 631, 904,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3590, 3, 652, 932,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3674, 3, 736,
                                                                       1016, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3758, 3, 757,
                                                                       1044, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3842, 3, 778,
                                                                       1072, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 3926, 3, 876,
                                                                       1172, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4034, 3, 904,
                                                                       1208, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4142, 3, 1016,
                                                                       1316, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 4250, 3, 1044,
                                                                       1352, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4358, 3, 1172,
                                                                       1478, ncols, p, q);

                    compute_prim_lsp_three_center_electron_repulsion_0(buffer, 4493, 3, 1316,
                                                                       1613, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4628, 3, 7, 8,
                                                                       1658, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4634, 3, 8, 9,
                                                                       1661, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4640, 3, 9, 10,
                                                                       1664, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4646, 3, 10, 11,
                                                                       1667, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4652, 3, 11, 12,
                                                                       1670, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4658, 3, 12, 13,
                                                                       1673, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4664, 3, 13, 14,
                                                                       1676, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4670, 3, 14, 15,
                                                                       1679, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4676, 3, 15, 16,
                                                                       1682, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4682, 3, 19, 20,
                                                                       1685, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4688, 3, 20, 21,
                                                                       1688, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4694, 3, 21, 22,
                                                                       1691, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4700, 3, 22, 23,
                                                                       1694, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4706, 3, 23, 24,
                                                                       1697, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4712, 3, 24, 25,
                                                                       1700, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4718, 3, 25, 26,
                                                                       1703, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4724, 3, 26, 27,
                                                                       1706, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4730, 3, 27, 28,
                                                                       1709, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4736, 0, 3, 4628,
                                                                       1658, 4634, 1712, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4754, 0, 3, 4634,
                                                                       1661, 4640, 1721, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4772, 0, 3, 4640,
                                                                       1664, 4646, 1730, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4790, 0, 3, 4646,
                                                                       1667, 4652, 1739, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4808, 0, 3, 4652,
                                                                       1670, 4658, 1748, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4826, 0, 3, 4658,
                                                                       1673, 4664, 1757, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4844, 0, 3, 4664,
                                                                       1676, 4670, 1766, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4862, 0, 3, 4670,
                                                                       1679, 4676, 1775, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4880, 0, 3, 4682,
                                                                       1685, 4688, 1784, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4898, 0, 3, 4688,
                                                                       1688, 4694, 1793, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4916, 0, 3, 4694,
                                                                       1691, 4700, 1802, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4934, 0, 3, 4700,
                                                                       1694, 4706, 1811, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4952, 0, 3, 4706,
                                                                       1697, 4712, 1820, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4970, 0, 3, 4712,
                                                                       1700, 4718, 1829, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4988, 0, 3, 4718,
                                                                       1703, 4724, 1838, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5006, 0, 3, 4724,
                                                                       1706, 4730, 1847, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5024, 0, 3, 4736,
                                                                       1712, 4754, 90, 96, 1856,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5060, 0, 3, 4754,
                                                                       1721, 4772, 96, 102, 1874,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5096, 0, 3, 4772,
                                                                       1730, 4790, 102, 108,
                                                                       1892, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5132, 0, 3, 4790,
                                                                       1739, 4808, 108, 114,
                                                                       1910, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5168, 0, 3, 4808,
                                                                       1748, 4826, 114, 120,
                                                                       1928, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5204, 0, 3, 4826,
                                                                       1757, 4844, 120, 126,
                                                                       1946, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5240, 0, 3, 4844,
                                                                       1766, 4862, 126, 132,
                                                                       1964, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5276, 0, 3, 4880,
                                                                       1784, 4898, 144, 150,
                                                                       1982, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5312, 0, 3, 4898,
                                                                       1793, 4916, 150, 156,
                                                                       2000, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5348, 0, 3, 4916,
                                                                       1802, 4934, 156, 162,
                                                                       2018, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5384, 0, 3, 4934,
                                                                       1811, 4952, 162, 168,
                                                                       2036, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5420, 0, 3, 4952,
                                                                       1820, 4970, 168, 174,
                                                                       2054, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5456, 0, 3, 4970,
                                                                       1829, 4988, 174, 180,
                                                                       2072, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5492, 0, 3, 4988,
                                                                       1838, 5006, 180, 186,
                                                                       2090, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5528, 0, 3, 5024,
                                                                       1856, 5060, 198, 208,
                                                                       2108, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5588, 0, 3, 5060,
                                                                       1874, 5096, 208, 218,
                                                                       2138, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5648, 0, 3, 5096,
                                                                       1892, 5132, 218, 228,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5708, 0, 3, 5132,
                                                                       1910, 5168, 228, 238,
                                                                       2198, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5768, 0, 3, 5168,
                                                                       1928, 5204, 238, 248,
                                                                       2228, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5828, 0, 3, 5204,
                                                                       1946, 5240, 248, 258,
                                                                       2258, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5888, 0, 3, 5276,
                                                                       1982, 5312, 278, 288,
                                                                       2288, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5948, 0, 3, 5312,
                                                                       2000, 5348, 288, 298,
                                                                       2318, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6008, 0, 3, 5348,
                                                                       2018, 5384, 298, 308,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6068, 0, 3, 5384,
                                                                       2036, 5420, 308, 318,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6128, 0, 3, 5420,
                                                                       2054, 5456, 318, 328,
                                                                       2408, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6188, 0, 3, 5456,
                                                                       2072, 5492, 328, 338,
                                                                       2438, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6248, 0, 3, 5528,
                                                                       2108, 5588, 358, 373,
                                                                       2468, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6338, 0, 3, 5588,
                                                                       2138, 5648, 373, 388,
                                                                       2513, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6428, 0, 3, 5648,
                                                                       2168, 5708, 388, 403,
                                                                       2558, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6518, 0, 3, 5708,
                                                                       2198, 5768, 403, 418,
                                                                       2603, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6608, 0, 3, 5768,
                                                                       2228, 5828, 418, 433,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6698, 0, 3, 5888,
                                                                       2288, 5948, 463, 478,
                                                                       2693, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6788, 0, 3, 5948,
                                                                       2318, 6008, 478, 493,
                                                                       2738, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6878, 0, 3, 6008,
                                                                       2348, 6068, 493, 508,
                                                                       2783, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6968, 0, 3, 6068,
                                                                       2378, 6128, 508, 523,
                                                                       2828, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7058, 0, 3, 6128,
                                                                       2408, 6188, 523, 538,
                                                                       2873, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7148, 0, 3, 6248,
                                                                       2468, 6338, 568, 589,
                                                                       2918, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7274, 0, 3, 6338,
                                                                       2513, 6428, 589, 610,
                                                                       2981, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 6428,
                                                                       2558, 6518, 610, 631,
                                                                       3044, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7526, 0, 3, 6518,
                                                                       2603, 6608, 631, 652,
                                                                       3107, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7652, 0, 3, 6698,
                                                                       2693, 6788, 694, 715,
                                                                       3170, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7778, 0, 3, 6788,
                                                                       2738, 6878, 715, 736,
                                                                       3233, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7904, 0, 3, 6878,
                                                                       2783, 6968, 736, 757,
                                                                       3296, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 8030, 0, 3, 6968,
                                                                       2828, 7058, 757, 778,
                                                                       3359, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8156, 0, 3, 7148,
                                                                       2918, 7274, 820, 848,
                                                                       3422, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8324, 0, 3, 7274,
                                                                       2981, 7400, 848, 876,
                                                                       3506, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8492, 0, 3, 7400,
                                                                       3044, 7526, 876, 904,
                                                                       3590, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8660, 0, 3, 7652,
                                                                       3170, 7778, 960, 988,
                                                                       3674, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8828, 0, 3, 7778,
                                                                       3233, 7904, 988, 1016,
                                                                       3758, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8996, 0, 3, 7904,
                                                                       3296, 8030, 1016, 1044,
                                                                       3842, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 9164, 0, 3, 8156,
                                                                       3422, 8324, 1100, 1136,
                                                                       3926, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 9380, 0, 3, 8324,
                                                                       3506, 8492, 1136, 1172,
                                                                       4034, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 9596, 0, 3, 8660,
                                                                       3674, 8828, 1244, 1280,
                                                                       4142, ncols, gamma, p,
                                                                       q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 9812, 0, 3, 8828,
                                                                       3758, 8996, 1280, 1316,
                                                                       4250, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 10028, 0, 3, 9164,
                                                                       3926, 9380, 1388, 1433,
                                                                       4358, ncols, gamma, p,
                                                                       q);

                    compute_prim_lsd_three_center_electron_repulsion_0(buffer, 10298, 0, 3, 9596,
                                                                       4142, 9812, 1523, 1568,
                                                                       4493, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_f_x(buffer, 10568, 5024, 6248, 1, 6, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 10628, 5024, 6248, 1, 6, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 10688, 5024, 6248, 1, 6, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 10748, 5276, 6698, 1, 6, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 10808, 5276, 6698, 1, 6, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 10868, 5276, 6698, 1, 6, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 10928, 5528, 7148, 1, 6, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 11018, 5528, 7148, 1, 6, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 11108, 5528, 7148, 1, 6, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 11198, 5888, 7652, 1, 6, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 11288, 5888, 7652, 1, 6, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 11378, 5888, 7652, 1, 6, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 11468, 6248, 8156, 1, 6, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 11594, 6248, 8156, 1, 6, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 11720, 6248, 8156, 1, 6, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 11846, 6698, 8660, 1, 6, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 11972, 6698, 8660, 1, 6, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 12098, 6698, 8660, 1, 6, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 12224, 7148, 9164, 1, 6, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 12392, 7148, 9164, 1, 6, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 12560, 7148, 9164, 1, 6, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 12728, 7652, 9596, 1, 6, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 12896, 7652, 9596, 1, 6, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 13064, 7652, 9596, 1, 6, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 13232, 8156, 10028, 1, 6, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 13448, 8156, 10028, 1, 6, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 13664, 8156, 10028, 1, 6, ncols, alpha);

                    simdgeo::geom_k_x(buffer, 13880, 8660, 10298, 1, 6, ncols, alpha);

                    simdgeo::geom_k_y(buffer, 14096, 8660, 10298, 1, 6, ncols, alpha);

                    simdgeo::geom_k_z(buffer, 14312, 8660, 10298, 1, 6, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 14528, 10568, 60, ncols);

                    simdfunc::contract_primitives(buffer, 14638, 10628, 60, ncols);

                    simdfunc::contract_primitives(buffer, 14748, 10688, 60, ncols);

                    simdfunc::contract_primitives(buffer, 14858, 5528, 60, ncols);

                    simdfunc::contract_primitives(buffer, 14968, 10748, 60, ncols);

                    simdfunc::contract_primitives(buffer, 15078, 10808, 60, ncols);

                    simdfunc::contract_primitives(buffer, 15188, 10868, 60, ncols);

                    simdfunc::contract_primitives(buffer, 15298, 5888, 60, ncols);

                    simdfunc::contract_primitives(buffer, 15408, 10928, 90, ncols);

                    simdfunc::contract_primitives(buffer, 15573, 11018, 90, ncols);

                    simdfunc::contract_primitives(buffer, 15738, 11108, 90, ncols);

                    simdfunc::contract_primitives(buffer, 15903, 6248, 90, ncols);

                    simdfunc::contract_primitives(buffer, 16068, 11198, 90, ncols);

                    simdfunc::contract_primitives(buffer, 16233, 11288, 90, ncols);

                    simdfunc::contract_primitives(buffer, 16398, 11378, 90, ncols);

                    simdfunc::contract_primitives(buffer, 16563, 6698, 90, ncols);

                    simdfunc::contract_primitives(buffer, 16728, 11468, 126, ncols);

                    simdfunc::contract_primitives(buffer, 16959, 11594, 126, ncols);

                    simdfunc::contract_primitives(buffer, 17190, 11720, 126, ncols);

                    simdfunc::contract_primitives(buffer, 17421, 7148, 126, ncols);

                    simdfunc::contract_primitives(buffer, 17652, 11846, 126, ncols);

                    simdfunc::contract_primitives(buffer, 17883, 11972, 126, ncols);

                    simdfunc::contract_primitives(buffer, 18114, 12098, 126, ncols);

                    simdfunc::contract_primitives(buffer, 18345, 7652, 126, ncols);

                    simdfunc::contract_primitives(buffer, 18576, 12224, 168, ncols);

                    simdfunc::contract_primitives(buffer, 18884, 12392, 168, ncols);

                    simdfunc::contract_primitives(buffer, 19192, 12560, 168, ncols);

                    simdfunc::contract_primitives(buffer, 19500, 8156, 168, ncols);

                    simdfunc::contract_primitives(buffer, 19808, 12728, 168, ncols);

                    simdfunc::contract_primitives(buffer, 20116, 12896, 168, ncols);

                    simdfunc::contract_primitives(buffer, 20424, 13064, 168, ncols);

                    simdfunc::contract_primitives(buffer, 20732, 8660, 168, ncols);

                    simdfunc::contract_primitives(buffer, 21040, 13232, 216, ncols);

                    simdfunc::contract_primitives(buffer, 21436, 13448, 216, ncols);

                    simdfunc::contract_primitives(buffer, 21832, 13664, 216, ncols);

                    simdfunc::contract_primitives(buffer, 22228, 13880, 216, ncols);

                    simdfunc::contract_primitives(buffer, 22624, 14096, 216, ncols);

                    simdfunc::contract_primitives(buffer, 23020, 14312, 216, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 14588, 14528, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 14698, 14638, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 14808, 14748, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 14918, 14858, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15028, 14968, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15138, 15078, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15248, 15188, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15358, 15298, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15498, 15408, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15663, 15573, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15828, 15738, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 15993, 15903, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 16158, 16068, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 16323, 16233, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 16488, 16398, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 16653, 16563, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 16854, 16728, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 17085, 16959, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 17316, 17190, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 17547, 17421, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 17778, 17652, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 18009, 17883, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 18240, 18114, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 18471, 18345, 21, 1, nmax);

        simdtrf::transform_d_inner(buffer, 18744, 18576, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 19052, 18884, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 19360, 19192, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 19668, 19500, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 19976, 19808, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 20284, 20116, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 20592, 20424, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 20900, 20732, 28, 1, nmax);

        simdtrf::transform_d_inner(buffer, 21256, 21040, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 21652, 21436, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 22048, 21832, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 22444, 22228, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 22840, 22624, 36, 1, nmax);

        simdtrf::transform_d_inner(buffer, 23236, 23020, 36, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 23416, 14588, 14918,
                                                       15498, 5, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 23566, 14698, 14918,
                                                       15663, 5, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 23716, 14808, 14918,
                                                       15828, 5, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 23866, 14918, 15993, 5, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 24016, 15028, 15358,
                                                       16158, 5, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 24166, 15138, 15358,
                                                       16323, 5, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 24316, 15248, 15358,
                                                       16488, 5, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 24466, 15358, 16653, 5, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 24616, 15498, 15993,
                                                       16854, 5, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 24841, 15663, 15993,
                                                       17085, 5, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 25066, 15828, 15993,
                                                       17316, 5, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 25291, 15993, 17547, 5, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 25516, 16158, 16653,
                                                       17778, 5, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 25741, 16323, 16653,
                                                       18009, 5, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 25966, 16488, 16653,
                                                       18240, 5, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 26191, 16653, 18471, 5, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 26416, 16854, 17547,
                                                       18744, 5, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 26731, 17085, 17547,
                                                       19052, 5, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 27046, 17316, 17547,
                                                       19360, 5, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 27361, 17547, 19668, 5, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 27676, 17778, 18471,
                                                       19976, 5, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 27991, 18009, 18471,
                                                       20284, 5, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 28306, 18240, 18471,
                                                       20592, 5, nmax);

        simdtrf::compute_hrr_hp_out_of_first(buffer, coordinates, 28621, 18471, 20900, 5, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 28936, 18744, 19668,
                                                       21256, 5, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 29356, 19052, 19668,
                                                       21652, 5, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 29776, 19360, 19668,
                                                       22048, 5, nmax);

        simdtrf::compute_hrr_geom_100x_ip_out_of_first(buffer, coordinates, 30196, 19976, 20900,
                                                       22444, 5, nmax);

        simdtrf::compute_hrr_geom_100y_ip_out_of_first(buffer, coordinates, 30616, 20284, 20900,
                                                       22840, 5, nmax);

        simdtrf::compute_hrr_geom_100z_ip_out_of_first(buffer, coordinates, 31036, 20592, 20900,
                                                       23236, 5, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 31456, 23416, 23866,
                                                       24616, 5, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 31756, 23566, 23866,
                                                       24841, 5, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 32056, 23716, 23866,
                                                       25066, 5, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 32356, 23866, 25291, 5, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 32656, 24016, 24466,
                                                       25516, 5, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 32956, 24166, 24466,
                                                       25741, 5, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 33256, 24316, 24466,
                                                       25966, 5, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 33556, 24466, 26191, 5, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 33856, 24616, 25291,
                                                       26416, 5, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 34306, 24841, 25291,
                                                       26731, 5, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 34756, 25066, 25291,
                                                       27046, 5, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 35206, 25291, 27361, 5, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 35656, 25516, 26191,
                                                       27676, 5, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 36106, 25741, 26191,
                                                       27991, 5, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 36556, 25966, 26191,
                                                       28306, 5, nmax);

        simdtrf::compute_hrr_gd_out_of_first(buffer, coordinates, 37006, 26191, 28621, 5, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 37456, 26416, 27361,
                                                       28936, 5, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 38086, 26731, 27361,
                                                       29356, 5, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 38716, 27046, 27361,
                                                       29776, 5, nmax);

        simdtrf::compute_hrr_geom_100x_hd_out_of_first(buffer, coordinates, 39346, 27676, 28621,
                                                       30196, 5, nmax);

        simdtrf::compute_hrr_geom_100y_hd_out_of_first(buffer, coordinates, 39976, 27991, 28621,
                                                       30616, 5, nmax);

        simdtrf::compute_hrr_geom_100z_hd_out_of_first(buffer, coordinates, 40606, 28306, 28621,
                                                       31036, 5, nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 41236, 31456, 32356,
                                                       33856, 5, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 41736, 31756, 32356,
                                                       34306, 5, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 42236, 32056, 32356,
                                                       34756, 5, nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 42736, 32356, 35206, 5, nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 43236, 32656, 33556,
                                                       35656, 5, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 43736, 32956, 33556,
                                                       36106, 5, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 44236, 33256, 33556,
                                                       36556, 5, nmax);

        simdtrf::compute_hrr_ff_out_of_first(buffer, coordinates, 44736, 33556, 37006, 5, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 45236, 33856, 35206,
                                                       37456, 5, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 45986, 34306, 35206,
                                                       38086, 5, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 46736, 34756, 35206,
                                                       38716, 5, nmax);

        simdtrf::compute_hrr_geom_100x_gf_out_of_first(buffer, coordinates, 47486, 35656, 37006,
                                                       39346, 5, nmax);

        simdtrf::compute_hrr_geom_100y_gf_out_of_first(buffer, coordinates, 48236, 36106, 37006,
                                                       39976, 5, nmax);

        simdtrf::compute_hrr_geom_100z_gf_out_of_first(buffer, coordinates, 48986, 36556, 37006,
                                                       40606, 5, nmax);

        simdtrf::compute_hrr_geom_100x_fg_out_of_first(buffer, coordinates, 49736, 41236, 42736,
                                                       45236, 5, nmax);

        simdtrf::compute_hrr_geom_100y_fg_out_of_first(buffer, coordinates, 50486, 41736, 42736,
                                                       45986, 5, nmax);

        simdtrf::compute_hrr_geom_100z_fg_out_of_first(buffer, coordinates, 51236, 42236, 42736,
                                                       46736, 5, nmax);

        simdtrf::compute_hrr_geom_100x_fg_out_of_first(buffer, coordinates, 51986, 43236, 44736,
                                                       47486, 5, nmax);

        simdtrf::compute_hrr_geom_100y_fg_out_of_first(buffer, coordinates, 52736, 43736, 44736,
                                                       48236, 5, nmax);

        simdtrf::compute_hrr_geom_100z_fg_out_of_first(buffer, coordinates, 53486, 44236, 44736,
                                                       48986, 5, nmax);

        simdtrf::transform_g_inner(buffer, 54236, 51986, 10, 5, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 54236, 45, nmax);

        simdtrf::transform_g_inner(buffer, 54236, 52736, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 54236,
                                   45, nmax);

        simdtrf::transform_g_inner(buffer, 54236, 53486, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 630 * nvalues + n * npairs, nvalues, buffer, 54236,
                                   45, nmax);

        simdtrf::transform_g_inner(buffer, 54236, 49736, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 945 * nvalues + n * npairs, nvalues, buffer, 54236,
                                   45, nmax);

        simdtrf::transform_g_inner(buffer, 54236, 50486, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 1260 * nvalues + n * npairs, nvalues, buffer, 54236,
                                   45, nmax);

        simdtrf::transform_g_inner(buffer, 54236, 51236, 10, 5, nmax);

        simdtrf::transform_f_outer(values + 1575 * nvalues + n * npairs, nvalues, buffer, 54236,
                                   45, nmax);
    }

    for (size_t m = 0; m < 1890; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
