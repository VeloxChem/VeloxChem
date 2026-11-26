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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecFDH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGeom100XFD.hpp"
#include "SimdTransferGeom100XFP.hpp"
#include "SimdTransferGeom100XGP.hpp"
#include "SimdTransferGeom100YFD.hpp"
#include "SimdTransferGeom100YFP.hpp"
#include "SimdTransferGeom100YGP.hpp"
#include "SimdTransferGeom100ZFD.hpp"
#include "SimdTransferGeom100ZFP.hpp"
#include "SimdTransferGeom100ZGP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_fdh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_fdh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 57940, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 2310 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 57940, 37388, 10201, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10, 11}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 18, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10, 11}, ncols, fj,
                                                        i * nprim_b + j, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1100, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1103, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1106, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1109, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1112, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1115, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1118, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1121, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1124, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1127, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1130, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1133, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1136, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1139, 3, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1142, 3, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1145, 3, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1148, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1151, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1154, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1157, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1160, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1163, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1166, 3, 7, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1175, 3, 8, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1184, 3, 9, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1193, 3, 10, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1202, 3, 11, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1211, 3, 12, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1220, 3, 13, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1229, 3, 14, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1238, 3, 15, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1247, 3, 16, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1256, 3, 19, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1265, 3, 20, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1274, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1283, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1292, 3, 23, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1301, 3, 24, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1310, 3, 25, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1319, 3, 26, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1328, 3, 27, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1337, 3, 28, 87,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1346, 3, 30, 90,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1364, 3, 33, 96,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1382, 3, 36, 102,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1400, 3, 39, 108,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1418, 3, 42, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1436, 3, 45, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1454, 3, 48, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1472, 3, 51, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1490, 3, 54, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1508, 3, 60, 144,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1526, 3, 63, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1544, 3, 66, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1562, 3, 69, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1580, 3, 72, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1598, 3, 75, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1616, 3, 78, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1634, 3, 81, 186,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1652, 3, 84, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1670, 3, 90, 198,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1700, 3, 96, 208,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1730, 3, 102, 218,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1760, 3, 108, 228,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1790, 3, 114, 238,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1820, 3, 120, 248,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1850, 3, 126, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1880, 3, 132, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1910, 3, 144, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1940, 3, 150, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1970, 3, 156, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2000, 3, 162, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2030, 3, 168, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2060, 3, 174, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2090, 3, 180, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2120, 3, 186, 348,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2150, 3, 198, 358,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2195, 3, 208, 373,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2240, 3, 218, 388,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2285, 3, 228, 403,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2330, 3, 238, 418,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2375, 3, 248, 433,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2420, 3, 258, 448,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2465, 3, 278, 463,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2510, 3, 288, 478,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2555, 3, 298, 493,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2600, 3, 308, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2645, 3, 318, 523,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2690, 3, 328, 538,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2735, 3, 338, 553,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2780, 3, 358, 568,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2843, 3, 373, 589,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2906, 3, 388, 610,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 2969, 3, 403, 631,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3032, 3, 418, 652,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3095, 3, 433, 673,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3158, 3, 463, 694,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3221, 3, 478, 715,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3284, 3, 493, 736,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3347, 3, 508, 757,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3410, 3, 523, 778,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3473, 3, 538, 799,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3536, 3, 568, 820,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3620, 3, 589, 848,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3704, 3, 610, 876,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3788, 3, 631, 904,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3872, 3, 652, 932,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 3956, 3, 694, 960,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4040, 3, 715, 988,
                                                                       ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4124, 3, 736,
                                                                       1016, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4208, 3, 757,
                                                                       1044, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4292, 3, 778,
                                                                       1072, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4376, 3, 7, 8,
                                                                       1106, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4382, 3, 8, 9,
                                                                       1109, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4388, 3, 9, 10,
                                                                       1112, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4394, 3, 10, 11,
                                                                       1115, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4400, 3, 11, 12,
                                                                       1118, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4406, 3, 12, 13,
                                                                       1121, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4412, 3, 13, 14,
                                                                       1124, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4418, 3, 14, 15,
                                                                       1127, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4424, 3, 15, 16,
                                                                       1130, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4430, 3, 19, 20,
                                                                       1139, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4436, 3, 20, 21,
                                                                       1142, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4442, 3, 21, 22,
                                                                       1145, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4448, 3, 22, 23,
                                                                       1148, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4454, 3, 23, 24,
                                                                       1151, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4460, 3, 24, 25,
                                                                       1154, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4466, 3, 25, 26,
                                                                       1157, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4472, 3, 26, 27,
                                                                       1160, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4478, 3, 27, 28,
                                                                       1163, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4484, 0, 3, 4376,
                                                                       1106, 4382, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4502, 0, 3, 4382,
                                                                       1109, 4388, 1193, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4520, 0, 3, 4388,
                                                                       1112, 4394, 1202, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4538, 0, 3, 4394,
                                                                       1115, 4400, 1211, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4556, 0, 3, 4400,
                                                                       1118, 4406, 1220, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4574, 0, 3, 4406,
                                                                       1121, 4412, 1229, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4592, 0, 3, 4412,
                                                                       1124, 4418, 1238, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4610, 0, 3, 4418,
                                                                       1127, 4424, 1247, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4628, 0, 3, 4430,
                                                                       1139, 4436, 1274, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4646, 0, 3, 4436,
                                                                       1142, 4442, 1283, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4664, 0, 3, 4442,
                                                                       1145, 4448, 1292, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4682, 0, 3, 4448,
                                                                       1148, 4454, 1301, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4700, 0, 3, 4454,
                                                                       1151, 4460, 1310, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4718, 0, 3, 4460,
                                                                       1154, 4466, 1319, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4736, 0, 3, 4466,
                                                                       1157, 4472, 1328, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4754, 0, 3, 4472,
                                                                       1160, 4478, 1337, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4772, 0, 3, 4484,
                                                                       1184, 4502, 90, 96, 1382,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4808, 0, 3, 4502,
                                                                       1193, 4520, 96, 102, 1400,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4844, 0, 3, 4520,
                                                                       1202, 4538, 102, 108,
                                                                       1418, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4880, 0, 3, 4538,
                                                                       1211, 4556, 108, 114,
                                                                       1436, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4916, 0, 3, 4556,
                                                                       1220, 4574, 114, 120,
                                                                       1454, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4952, 0, 3, 4574,
                                                                       1229, 4592, 120, 126,
                                                                       1472, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4988, 0, 3, 4592,
                                                                       1238, 4610, 126, 132,
                                                                       1490, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5024, 0, 3, 4628,
                                                                       1274, 4646, 144, 150,
                                                                       1544, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5060, 0, 3, 4646,
                                                                       1283, 4664, 150, 156,
                                                                       1562, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5096, 0, 3, 4664,
                                                                       1292, 4682, 156, 162,
                                                                       1580, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5132, 0, 3, 4682,
                                                                       1301, 4700, 162, 168,
                                                                       1598, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5168, 0, 3, 4700,
                                                                       1310, 4718, 168, 174,
                                                                       1616, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5204, 0, 3, 4718,
                                                                       1319, 4736, 174, 180,
                                                                       1634, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5240, 0, 3, 4736,
                                                                       1328, 4754, 180, 186,
                                                                       1652, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5276, 0, 3, 4772,
                                                                       1382, 4808, 198, 208,
                                                                       1730, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5336, 0, 3, 4808,
                                                                       1400, 4844, 208, 218,
                                                                       1760, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5396, 0, 3, 4844,
                                                                       1418, 4880, 218, 228,
                                                                       1790, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5456, 0, 3, 4880,
                                                                       1436, 4916, 228, 238,
                                                                       1820, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5516, 0, 3, 4916,
                                                                       1454, 4952, 238, 248,
                                                                       1850, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5576, 0, 3, 4952,
                                                                       1472, 4988, 248, 258,
                                                                       1880, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5636, 0, 3, 5024,
                                                                       1544, 5060, 278, 288,
                                                                       1970, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5696, 0, 3, 5060,
                                                                       1562, 5096, 288, 298,
                                                                       2000, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5756, 0, 3, 5096,
                                                                       1580, 5132, 298, 308,
                                                                       2030, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5816, 0, 3, 5132,
                                                                       1598, 5168, 308, 318,
                                                                       2060, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5876, 0, 3, 5168,
                                                                       1616, 5204, 318, 328,
                                                                       2090, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5936, 0, 3, 5204,
                                                                       1634, 5240, 328, 338,
                                                                       2120, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5996, 0, 3, 5276,
                                                                       1730, 5336, 358, 373,
                                                                       2240, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6086, 0, 3, 5336,
                                                                       1760, 5396, 373, 388,
                                                                       2285, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6176, 0, 3, 5396,
                                                                       1790, 5456, 388, 403,
                                                                       2330, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6266, 0, 3, 5456,
                                                                       1820, 5516, 403, 418,
                                                                       2375, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6356, 0, 3, 5516,
                                                                       1850, 5576, 418, 433,
                                                                       2420, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6446, 0, 3, 5636,
                                                                       1970, 5696, 463, 478,
                                                                       2555, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6536, 0, 3, 5696,
                                                                       2000, 5756, 478, 493,
                                                                       2600, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6626, 0, 3, 5756,
                                                                       2030, 5816, 493, 508,
                                                                       2645, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6716, 0, 3, 5816,
                                                                       2060, 5876, 508, 523,
                                                                       2690, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6806, 0, 3, 5876,
                                                                       2090, 5936, 523, 538,
                                                                       2735, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 6896, 0, 3, 5996,
                                                                       2240, 6086, 568, 589,
                                                                       2906, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7022, 0, 3, 6086,
                                                                       2285, 6176, 589, 610,
                                                                       2969, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7148, 0, 3, 6176,
                                                                       2330, 6266, 610, 631,
                                                                       3032, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7274, 0, 3, 6266,
                                                                       2375, 6356, 631, 652,
                                                                       3095, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 6446,
                                                                       2555, 6536, 694, 715,
                                                                       3284, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7526, 0, 3, 6536,
                                                                       2600, 6626, 715, 736,
                                                                       3347, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7652, 0, 3, 6626,
                                                                       2645, 6716, 736, 757,
                                                                       3410, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 7778, 0, 3, 6716,
                                                                       2690, 6806, 757, 778,
                                                                       3473, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 7904, 0, 3, 6896,
                                                                       2906, 7022, 820, 848,
                                                                       3704, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8072, 0, 3, 7022,
                                                                       2969, 7148, 848, 876,
                                                                       3788, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8240, 0, 3, 7148,
                                                                       3032, 7274, 876, 904,
                                                                       3872, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8408, 0, 3, 7400,
                                                                       3284, 7526, 960, 988,
                                                                       4124, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8576, 0, 3, 7526,
                                                                       3347, 7652, 988, 1016,
                                                                       4208, ncols, gamma, p,
                                                                       q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 8744, 0, 3, 7652,
                                                                       3410, 7778, 1016, 1044,
                                                                       4292, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8912, 3, 1100,
                                                                       1103, 4376, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8922, 3, 1103,
                                                                       1106, 4382, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8932, 3, 1106,
                                                                       1109, 4388, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8942, 3, 1109,
                                                                       1112, 4394, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8952, 3, 1112,
                                                                       1115, 4400, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8962, 3, 1115,
                                                                       1118, 4406, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8972, 3, 1118,
                                                                       1121, 4412, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8982, 3, 1121,
                                                                       1124, 4418, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8992, 3, 1124,
                                                                       1127, 4424, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9002, 3, 1133,
                                                                       1136, 4430, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9012, 3, 1136,
                                                                       1139, 4436, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9022, 3, 1139,
                                                                       1142, 4442, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9032, 3, 1142,
                                                                       1145, 4448, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9042, 3, 1145,
                                                                       1148, 4454, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9052, 3, 1148,
                                                                       1151, 4460, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9062, 3, 1151,
                                                                       1154, 4466, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9072, 3, 1154,
                                                                       1157, 4472, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9082, 3, 1157,
                                                                       1160, 4478, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9092, 0, 3, 8912,
                                                                       4376, 8922, 1166, 1175,
                                                                       4484, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9122, 0, 3, 8922,
                                                                       4382, 8932, 1175, 1184,
                                                                       4502, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9152, 0, 3, 8932,
                                                                       4388, 8942, 1184, 1193,
                                                                       4520, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9182, 0, 3, 8942,
                                                                       4394, 8952, 1193, 1202,
                                                                       4538, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9212, 0, 3, 8952,
                                                                       4400, 8962, 1202, 1211,
                                                                       4556, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9242, 0, 3, 8962,
                                                                       4406, 8972, 1211, 1220,
                                                                       4574, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9272, 0, 3, 8972,
                                                                       4412, 8982, 1220, 1229,
                                                                       4592, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9302, 0, 3, 8982,
                                                                       4418, 8992, 1229, 1238,
                                                                       4610, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9332, 0, 3, 9002,
                                                                       4430, 9012, 1256, 1265,
                                                                       4628, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9362, 0, 3, 9012,
                                                                       4436, 9022, 1265, 1274,
                                                                       4646, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9392, 0, 3, 9022,
                                                                       4442, 9032, 1274, 1283,
                                                                       4664, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9422, 0, 3, 9032,
                                                                       4448, 9042, 1283, 1292,
                                                                       4682, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9452, 0, 3, 9042,
                                                                       4454, 9052, 1292, 1301,
                                                                       4700, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9482, 0, 3, 9052,
                                                                       4460, 9062, 1301, 1310,
                                                                       4718, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9512, 0, 3, 9062,
                                                                       4466, 9072, 1310, 1319,
                                                                       4736, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9542, 0, 3, 9072,
                                                                       4472, 9082, 1319, 1328,
                                                                       4754, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9572, 0, 3, 9092,
                                                                       4484, 9122, 1346, 1364,
                                                                       4772, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9632, 0, 3, 9122,
                                                                       4502, 9152, 1364, 1382,
                                                                       4808, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9692, 0, 3, 9152,
                                                                       4520, 9182, 1382, 1400,
                                                                       4844, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9752, 0, 3, 9182,
                                                                       4538, 9212, 1400, 1418,
                                                                       4880, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9812, 0, 3, 9212,
                                                                       4556, 9242, 1418, 1436,
                                                                       4916, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9872, 0, 3, 9242,
                                                                       4574, 9272, 1436, 1454,
                                                                       4952, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9932, 0, 3, 9272,
                                                                       4592, 9302, 1454, 1472,
                                                                       4988, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9992, 0, 3, 9332,
                                                                       4628, 9362, 1508, 1526,
                                                                       5024, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10052, 0, 3, 9362,
                                                                       4646, 9392, 1526, 1544,
                                                                       5060, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10112, 0, 3, 9392,
                                                                       4664, 9422, 1544, 1562,
                                                                       5096, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10172, 0, 3, 9422,
                                                                       4682, 9452, 1562, 1580,
                                                                       5132, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10232, 0, 3, 9452,
                                                                       4700, 9482, 1580, 1598,
                                                                       5168, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10292, 0, 3, 9482,
                                                                       4718, 9512, 1598, 1616,
                                                                       5204, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10352, 0, 3, 9512,
                                                                       4736, 9542, 1616, 1634,
                                                                       5240, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10412, 0, 3, 9572,
                                                                       4772, 9632, 1670, 1700,
                                                                       5276, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10512, 0, 3, 9632,
                                                                       4808, 9692, 1700, 1730,
                                                                       5336, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10612, 0, 3, 9692,
                                                                       4844, 9752, 1730, 1760,
                                                                       5396, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10712, 0, 3, 9752,
                                                                       4880, 9812, 1760, 1790,
                                                                       5456, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10812, 0, 3, 9812,
                                                                       4916, 9872, 1790, 1820,
                                                                       5516, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 10912, 0, 3, 9872,
                                                                       4952, 9932, 1820, 1850,
                                                                       5576, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11012, 0, 3, 9992,
                                                                       5024, 10052, 1910, 1940,
                                                                       5636, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11112, 0, 3,
                                                                       10052, 5060, 10112, 1940,
                                                                       1970, 5696, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11212, 0, 3,
                                                                       10112, 5096, 10172, 1970,
                                                                       2000, 5756, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11312, 0, 3,
                                                                       10172, 5132, 10232, 2000,
                                                                       2030, 5816, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11412, 0, 3,
                                                                       10232, 5168, 10292, 2030,
                                                                       2060, 5876, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11512, 0, 3,
                                                                       10292, 5204, 10352, 2060,
                                                                       2090, 5936, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11612, 0, 3,
                                                                       10412, 5276, 10512, 2150,
                                                                       2195, 5996, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11762, 0, 3,
                                                                       10512, 5336, 10612, 2195,
                                                                       2240, 6086, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 11912, 0, 3,
                                                                       10612, 5396, 10712, 2240,
                                                                       2285, 6176, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 12062, 0, 3,
                                                                       10712, 5456, 10812, 2285,
                                                                       2330, 6266, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 12212, 0, 3,
                                                                       10812, 5516, 10912, 2330,
                                                                       2375, 6356, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 12362, 0, 3,
                                                                       11012, 5636, 11112, 2465,
                                                                       2510, 6446, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 12512, 0, 3,
                                                                       11112, 5696, 11212, 2510,
                                                                       2555, 6536, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 12662, 0, 3,
                                                                       11212, 5756, 11312, 2555,
                                                                       2600, 6626, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 12812, 0, 3,
                                                                       11312, 5816, 11412, 2600,
                                                                       2645, 6716, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 12962, 0, 3,
                                                                       11412, 5876, 11512, 2645,
                                                                       2690, 6806, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 13112, 0, 3,
                                                                       11612, 5996, 11762, 2780,
                                                                       2843, 6896, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 13322, 0, 3,
                                                                       11762, 6086, 11912, 2843,
                                                                       2906, 7022, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 13532, 0, 3,
                                                                       11912, 6176, 12062, 2906,
                                                                       2969, 7148, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 13742, 0, 3,
                                                                       12062, 6266, 12212, 2969,
                                                                       3032, 7274, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 13952, 0, 3,
                                                                       12362, 6446, 12512, 3158,
                                                                       3221, 7400, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 14162, 0, 3,
                                                                       12512, 6536, 12662, 3221,
                                                                       3284, 7526, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 14372, 0, 3,
                                                                       12662, 6626, 12812, 3284,
                                                                       3347, 7652, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 14582, 0, 3,
                                                                       12812, 6716, 12962, 3347,
                                                                       3410, 7778, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 14792, 0, 3,
                                                                       13112, 6896, 13322, 3536,
                                                                       3620, 7904, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 15072, 0, 3,
                                                                       13322, 7022, 13532, 3620,
                                                                       3704, 8072, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 15352, 0, 3,
                                                                       13532, 7148, 13742, 3704,
                                                                       3788, 8240, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 15632, 0, 3,
                                                                       13952, 7400, 14162, 3956,
                                                                       4040, 8408, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 15912, 0, 3,
                                                                       14162, 7526, 14372, 4040,
                                                                       4124, 8576, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 16192, 0, 3,
                                                                       14372, 7652, 14582, 4124,
                                                                       4208, 8744, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16472, 3, 4376,
                                                                       4382, 8932, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16487, 3, 4382,
                                                                       4388, 8942, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16502, 3, 4388,
                                                                       4394, 8952, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16517, 3, 4394,
                                                                       4400, 8962, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16532, 3, 4400,
                                                                       4406, 8972, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16547, 3, 4406,
                                                                       4412, 8982, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16562, 3, 4412,
                                                                       4418, 8992, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16577, 3, 4430,
                                                                       4436, 9022, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16592, 3, 4436,
                                                                       4442, 9032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16607, 3, 4442,
                                                                       4448, 9042, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16622, 3, 4448,
                                                                       4454, 9052, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16637, 3, 4454,
                                                                       4460, 9062, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16652, 3, 4460,
                                                                       4466, 9072, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16667, 3, 4466,
                                                                       4472, 9082, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 16682, 0, 3,
                                                                       16472, 8932, 16487, 4484,
                                                                       4502, 9152, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 16727, 0, 3,
                                                                       16487, 8942, 16502, 4502,
                                                                       4520, 9182, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 16772, 0, 3,
                                                                       16502, 8952, 16517, 4520,
                                                                       4538, 9212, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 16817, 0, 3,
                                                                       16517, 8962, 16532, 4538,
                                                                       4556, 9242, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 16862, 0, 3,
                                                                       16532, 8972, 16547, 4556,
                                                                       4574, 9272, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 16907, 0, 3,
                                                                       16547, 8982, 16562, 4574,
                                                                       4592, 9302, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 16952, 0, 3,
                                                                       16577, 9022, 16592, 4628,
                                                                       4646, 9392, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 16997, 0, 3,
                                                                       16592, 9032, 16607, 4646,
                                                                       4664, 9422, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17042, 0, 3,
                                                                       16607, 9042, 16622, 4664,
                                                                       4682, 9452, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17087, 0, 3,
                                                                       16622, 9052, 16637, 4682,
                                                                       4700, 9482, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17132, 0, 3,
                                                                       16637, 9062, 16652, 4700,
                                                                       4718, 9512, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17177, 0, 3,
                                                                       16652, 9072, 16667, 4718,
                                                                       4736, 9542, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17222, 0, 3,
                                                                       16682, 9152, 16727, 4772,
                                                                       4808, 9692, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17312, 0, 3,
                                                                       16727, 9182, 16772, 4808,
                                                                       4844, 9752, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17402, 0, 3,
                                                                       16772, 9212, 16817, 4844,
                                                                       4880, 9812, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17492, 0, 3,
                                                                       16817, 9242, 16862, 4880,
                                                                       4916, 9872, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17582, 0, 3,
                                                                       16862, 9272, 16907, 4916,
                                                                       4952, 9932, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17672, 0, 3,
                                                                       16952, 9392, 16997, 5024,
                                                                       5060, 10112, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17762, 0, 3,
                                                                       16997, 9422, 17042, 5060,
                                                                       5096, 10172, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17852, 0, 3,
                                                                       17042, 9452, 17087, 5096,
                                                                       5132, 10232, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17942, 0, 3,
                                                                       17087, 9482, 17132, 5132,
                                                                       5168, 10292, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 18032, 0, 3,
                                                                       17132, 9512, 17177, 5168,
                                                                       5204, 10352, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18122, 0, 3,
                                                                       17222, 9692, 17312, 5276,
                                                                       5336, 10612, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18272, 0, 3,
                                                                       17312, 9752, 17402, 5336,
                                                                       5396, 10712, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18422, 0, 3,
                                                                       17402, 9812, 17492, 5396,
                                                                       5456, 10812, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18572, 0, 3,
                                                                       17492, 9872, 17582, 5456,
                                                                       5516, 10912, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18722, 0, 3,
                                                                       17672, 10112, 17762, 5636,
                                                                       5696, 11212, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18872, 0, 3,
                                                                       17762, 10172, 17852, 5696,
                                                                       5756, 11312, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 19022, 0, 3,
                                                                       17852, 10232, 17942, 5756,
                                                                       5816, 11412, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 19172, 0, 3,
                                                                       17942, 10292, 18032, 5816,
                                                                       5876, 11512, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 19322, 0, 3,
                                                                       18122, 10612, 18272, 5996,
                                                                       6086, 11912, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 19547, 0, 3,
                                                                       18272, 10712, 18422, 6086,
                                                                       6176, 12062, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 19772, 0, 3,
                                                                       18422, 10812, 18572, 6176,
                                                                       6266, 12212, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 19997, 0, 3,
                                                                       18722, 11212, 18872, 6446,
                                                                       6536, 12662, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 20222, 0, 3,
                                                                       18872, 11312, 19022, 6536,
                                                                       6626, 12812, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 20447, 0, 3,
                                                                       19022, 11412, 19172, 6626,
                                                                       6716, 12962, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 20672, 0, 3,
                                                                       19322, 11912, 19547, 6896,
                                                                       7022, 13532, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 20987, 0, 3,
                                                                       19547, 12062, 19772, 7022,
                                                                       7148, 13742, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 21302, 0, 3,
                                                                       19997, 12662, 20222, 7400,
                                                                       7526, 14372, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 21617, 0, 3,
                                                                       20222, 12812, 20447, 7526,
                                                                       7652, 14582, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 21932, 0, 3,
                                                                       20672, 13532, 20987, 7904,
                                                                       8072, 15352, ncols, gamma,
                                                                       p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 22352, 0, 3,
                                                                       21302, 14372, 21617, 8408,
                                                                       8576, 16192, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22772, 3, 8912,
                                                                       8922, 16472, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22793, 3, 8922,
                                                                       8932, 16487, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22814, 3, 8932,
                                                                       8942, 16502, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22835, 3, 8942,
                                                                       8952, 16517, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22856, 3, 8952,
                                                                       8962, 16532, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22877, 3, 8962,
                                                                       8972, 16547, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22898, 3, 8972,
                                                                       8982, 16562, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22919, 3, 9002,
                                                                       9012, 16577, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22940, 3, 9012,
                                                                       9022, 16592, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22961, 3, 9022,
                                                                       9032, 16607, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 22982, 3, 9032,
                                                                       9042, 16622, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 23003, 3, 9042,
                                                                       9052, 16637, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 23024, 3, 9052,
                                                                       9062, 16652, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 23045, 3, 9062,
                                                                       9072, 16667, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23066, 0, 3,
                                                                       22772, 16472, 22793, 9092,
                                                                       9122, 16682, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23129, 0, 3,
                                                                       22793, 16487, 22814, 9122,
                                                                       9152, 16727, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23192, 0, 3,
                                                                       22814, 16502, 22835, 9152,
                                                                       9182, 16772, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23255, 0, 3,
                                                                       22835, 16517, 22856, 9182,
                                                                       9212, 16817, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23318, 0, 3,
                                                                       22856, 16532, 22877, 9212,
                                                                       9242, 16862, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23381, 0, 3,
                                                                       22877, 16547, 22898, 9242,
                                                                       9272, 16907, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23444, 0, 3,
                                                                       22919, 16577, 22940, 9332,
                                                                       9362, 16952, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23507, 0, 3,
                                                                       22940, 16592, 22961, 9362,
                                                                       9392, 16997, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23570, 0, 3,
                                                                       22961, 16607, 22982, 9392,
                                                                       9422, 17042, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23633, 0, 3,
                                                                       22982, 16622, 23003, 9422,
                                                                       9452, 17087, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23696, 0, 3,
                                                                       23003, 16637, 23024, 9452,
                                                                       9482, 17132, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 23759, 0, 3,
                                                                       23024, 16652, 23045, 9482,
                                                                       9512, 17177, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23822, 0, 3,
                                                                       23066, 16682, 23129, 9572,
                                                                       9632, 17222, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 23948, 0, 3,
                                                                       23129, 16727, 23192, 9632,
                                                                       9692, 17312, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 24074, 0, 3,
                                                                       23192, 16772, 23255, 9692,
                                                                       9752, 17402, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 24200, 0, 3,
                                                                       23255, 16817, 23318, 9752,
                                                                       9812, 17492, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 24326, 0, 3,
                                                                       23318, 16862, 23381, 9812,
                                                                       9872, 17582, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 24452, 0, 3,
                                                                       23444, 16952, 23507, 9992,
                                                                       10052, 17672, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 24578, 0, 3,
                                                                       23507, 16997, 23570,
                                                                       10052, 10112, 17762,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 24704, 0, 3,
                                                                       23570, 17042, 23633,
                                                                       10112, 10172, 17852,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 24830, 0, 3,
                                                                       23633, 17087, 23696,
                                                                       10172, 10232, 17942,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 24956, 0, 3,
                                                                       23696, 17132, 23759,
                                                                       10232, 10292, 18032,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 25082, 0, 3,
                                                                       23822, 17222, 23948,
                                                                       10412, 10512, 18122,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 25292, 0, 3,
                                                                       23948, 17312, 24074,
                                                                       10512, 10612, 18272,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 25502, 0, 3,
                                                                       24074, 17402, 24200,
                                                                       10612, 10712, 18422,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 25712, 0, 3,
                                                                       24200, 17492, 24326,
                                                                       10712, 10812, 18572,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 25922, 0, 3,
                                                                       24452, 17672, 24578,
                                                                       11012, 11112, 18722,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 26132, 0, 3,
                                                                       24578, 17762, 24704,
                                                                       11112, 11212, 18872,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 26342, 0, 3,
                                                                       24704, 17852, 24830,
                                                                       11212, 11312, 19022,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 26552, 0, 3,
                                                                       24830, 17942, 24956,
                                                                       11312, 11412, 19172,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 26762, 0, 3,
                                                                       25082, 18122, 25292,
                                                                       11612, 11762, 19322,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 27077, 0, 3,
                                                                       25292, 18272, 25502,
                                                                       11762, 11912, 19547,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 27392, 0, 3,
                                                                       25502, 18422, 25712,
                                                                       11912, 12062, 19772,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 27707, 0, 3,
                                                                       25922, 18722, 26132,
                                                                       12362, 12512, 19997,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 28022, 0, 3,
                                                                       26132, 18872, 26342,
                                                                       12512, 12662, 20222,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 28337, 0, 3,
                                                                       26342, 19022, 26552,
                                                                       12662, 12812, 20447,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 28652, 0, 3,
                                                                       26762, 19322, 27077,
                                                                       13112, 13322, 20672,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 29093, 0, 3,
                                                                       27077, 19547, 27392,
                                                                       13322, 13532, 20987,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 29534, 0, 3,
                                                                       27707, 19997, 28022,
                                                                       13952, 14162, 21302,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 29975, 0, 3,
                                                                       28022, 20222, 28337,
                                                                       14162, 14372, 21617,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 30416, 0, 3,
                                                                       28652, 20672, 29093,
                                                                       14792, 15072, 21932,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 31004, 0, 3,
                                                                       29534, 21302, 29975,
                                                                       15632, 15912, 22352,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_f_x(buffer, 31592, 23822, 26762, 1, 21, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 31802, 23822, 26762, 1, 21, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 32012, 23822, 26762, 1, 21, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 32222, 24452, 27707, 1, 21, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 32432, 24452, 27707, 1, 21, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 32642, 24452, 27707, 1, 21, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 32852, 25082, 28652, 1, 21, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 33167, 25082, 28652, 1, 21, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 33482, 25082, 28652, 1, 21, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 33797, 25922, 29534, 1, 21, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 34112, 25922, 29534, 1, 21, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 34427, 25922, 29534, 1, 21, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 34742, 26762, 30416, 1, 21, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 35183, 26762, 30416, 1, 21, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 35624, 26762, 30416, 1, 21, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 36065, 27707, 31004, 1, 21, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 36506, 27707, 31004, 1, 21, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 36947, 27707, 31004, 1, 21, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 37388, 31592, 210, ncols);

                    simdfunc::contract_primitives(buffer, 37708, 31802, 210, ncols);

                    simdfunc::contract_primitives(buffer, 38028, 32012, 210, ncols);

                    simdfunc::contract_primitives(buffer, 38348, 25082, 210, ncols);

                    simdfunc::contract_primitives(buffer, 38668, 32222, 210, ncols);

                    simdfunc::contract_primitives(buffer, 38988, 32432, 210, ncols);

                    simdfunc::contract_primitives(buffer, 39308, 32642, 210, ncols);

                    simdfunc::contract_primitives(buffer, 39628, 25922, 210, ncols);

                    simdfunc::contract_primitives(buffer, 39948, 32852, 315, ncols);

                    simdfunc::contract_primitives(buffer, 40428, 33167, 315, ncols);

                    simdfunc::contract_primitives(buffer, 40908, 33482, 315, ncols);

                    simdfunc::contract_primitives(buffer, 41388, 26762, 315, ncols);

                    simdfunc::contract_primitives(buffer, 41868, 33797, 315, ncols);

                    simdfunc::contract_primitives(buffer, 42348, 34112, 315, ncols);

                    simdfunc::contract_primitives(buffer, 42828, 34427, 315, ncols);

                    simdfunc::contract_primitives(buffer, 43308, 27707, 315, ncols);

                    simdfunc::contract_primitives(buffer, 43788, 34742, 441, ncols);

                    simdfunc::contract_primitives(buffer, 44460, 35183, 441, ncols);

                    simdfunc::contract_primitives(buffer, 45132, 35624, 441, ncols);

                    simdfunc::contract_primitives(buffer, 45804, 36065, 441, ncols);

                    simdfunc::contract_primitives(buffer, 46476, 36506, 441, ncols);

                    simdfunc::contract_primitives(buffer, 47148, 36947, 441, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 37598, 37388, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 37918, 37708, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 38238, 38028, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 38558, 38348, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 38878, 38668, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 39198, 38988, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 39518, 39308, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 39838, 39628, 10, 1, nmax);

        simdtrf::transform_h_inner(buffer, 40263, 39948, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 40743, 40428, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 41223, 40908, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 41703, 41388, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 42183, 41868, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 42663, 42348, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 43143, 42828, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 43623, 43308, 15, 1, nmax);

        simdtrf::transform_h_inner(buffer, 44229, 43788, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 44901, 44460, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 45573, 45132, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 46245, 45804, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 46917, 46476, 21, 1, nmax);

        simdtrf::transform_h_inner(buffer, 47589, 47148, 21, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 47820, 37598, 38558,
                                                       40263, 11, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 48150, 37918, 38558,
                                                       40743, 11, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 48480, 38238, 38558,
                                                       41223, 11, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 48810, 38558, 41703, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 49140, 38878, 39838,
                                                       42183, 11, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 49470, 39198, 39838,
                                                       42663, 11, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 49800, 39518, 39838,
                                                       43143, 11, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 50130, 39838, 43623, 11,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 50460, 40263, 41703,
                                                       44229, 11, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 50955, 40743, 41703,
                                                       44901, 11, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 51450, 41223, 41703,
                                                       45573, 11, nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 51945, 42183, 43623,
                                                       46245, 11, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 52440, 42663, 43623,
                                                       46917, 11, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 52935, 43143, 43623,
                                                       47589, 11, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 53430, 47820, 48810,
                                                       50460, 11, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 54090, 48150, 48810,
                                                       50955, 11, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 54750, 48480, 48810,
                                                       51450, 11, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 55410, 49140, 50130,
                                                       51945, 11, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 56070, 49470, 50130,
                                                       52440, 11, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 56730, 49800, 50130,
                                                       52935, 11, nmax);

        simdtrf::transform_d_inner(buffer, 57390, 55410, 10, 11, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 57390, 55, nmax);

        simdtrf::transform_d_inner(buffer, 57390, 56070, 10, 11, nmax);

        simdtrf::transform_f_outer(values + 385 * nvalues + n * npairs, nvalues, buffer, 57390,
                                   55, nmax);

        simdtrf::transform_d_inner(buffer, 57390, 56730, 10, 11, nmax);

        simdtrf::transform_f_outer(values + 770 * nvalues + n * npairs, nvalues, buffer, 57390,
                                   55, nmax);

        simdtrf::transform_d_inner(buffer, 57390, 53430, 10, 11, nmax);

        simdtrf::transform_f_outer(values + 1155 * nvalues + n * npairs, nvalues, buffer, 57390,
                                   55, nmax);

        simdtrf::transform_d_inner(buffer, 57390, 54090, 10, 11, nmax);

        simdtrf::transform_f_outer(values + 1540 * nvalues + n * npairs, nvalues, buffer, 57390,
                                   55, nmax);

        simdtrf::transform_d_inner(buffer, 57390, 54750, 10, 11, nmax);

        simdtrf::transform_f_outer(values + 1925 * nvalues + n * npairs, nvalues, buffer, 57390,
                                   55, nmax);
    }

    for (size_t m = 0; m < 2310; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
