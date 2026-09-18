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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecFFI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecISS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecKSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferFD.hpp"
#include "SimdTransferFP.hpp"
#include "SimdTransferGP.hpp"
#include "SimdTransferGeom100XFD.hpp"
#include "SimdTransferGeom100XFF.hpp"
#include "SimdTransferGeom100XFP.hpp"
#include "SimdTransferGeom100XGD.hpp"
#include "SimdTransferGeom100XGP.hpp"
#include "SimdTransferGeom100XHP.hpp"
#include "SimdTransferGeom100YFD.hpp"
#include "SimdTransferGeom100YFF.hpp"
#include "SimdTransferGeom100YFP.hpp"
#include "SimdTransferGeom100YGD.hpp"
#include "SimdTransferGeom100YGP.hpp"
#include "SimdTransferGeom100YHP.hpp"
#include "SimdTransferGeom100ZFD.hpp"
#include "SimdTransferGeom100ZFF.hpp"
#include "SimdTransferGeom100ZFP.hpp"
#include "SimdTransferGeom100ZGD.hpp"
#include "SimdTransferGeom100ZGP.hpp"
#include "SimdTransferGeom100ZHP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_ffi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_ffi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 146620, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 3822 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 146620, 89960, 21612, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 13,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 21, 3, 13,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 60, 0, 3, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 63, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 66, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 69, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 75, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 78, 0, 3, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 81, 0, 3, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 84, 0, 3, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 87, 0, 3, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 90, 0, 3, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 93, 0, 3, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 96, 0, 3, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 99, 0, 3, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 105, 0, 3, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 111, 0, 3, 34, 35,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 114, 0, 3, 7, 8,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 120, 0, 3, 8, 9,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 126, 0, 3, 9, 10,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 10, 11,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 11, 12,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 144, 0, 3, 12, 13,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 150, 0, 3, 13, 14,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 156, 0, 3, 14, 15,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 162, 0, 3, 15, 16,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 168, 0, 3, 16, 17,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 174, 0, 3, 17, 18,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 180, 0, 3, 18, 19,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 186, 0, 3, 22, 23,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 192, 0, 3, 23, 24,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 198, 0, 3, 24, 25,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 204, 0, 3, 25, 26,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 210, 0, 3, 26, 27,
                                                                       87, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 216, 0, 3, 27, 28,
                                                                       90, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 222, 0, 3, 28, 29,
                                                                       93, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 228, 0, 3, 29, 30,
                                                                       96, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 234, 0, 3, 30, 31,
                                                                       99, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 240, 0, 3, 31, 32,
                                                                       102, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 246, 0, 3, 32, 33,
                                                                       105, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 252, 0, 3, 33, 34,
                                                                       108, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 258, 0, 3, 36, 39,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 39, 42,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 278, 0, 3, 42, 45,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 288, 0, 3, 45, 48,
                                                                       132, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 48, 51,
                                                                       138, 144, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 308, 0, 3, 51, 54,
                                                                       144, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 318, 0, 3, 54, 57,
                                                                       150, 156, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 57, 60,
                                                                       156, 162, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 338, 0, 3, 60, 63,
                                                                       162, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 348, 0, 3, 63, 66,
                                                                       168, 174, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 66, 69,
                                                                       174, 180, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 368, 0, 3, 75, 78,
                                                                       186, 192, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 378, 0, 3, 78, 81,
                                                                       192, 198, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 388, 0, 3, 81, 84,
                                                                       198, 204, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 398, 0, 3, 84, 87,
                                                                       204, 210, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 408, 0, 3, 87, 90,
                                                                       210, 216, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 418, 0, 3, 90, 93,
                                                                       216, 222, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 428, 0, 3, 93, 96,
                                                                       222, 228, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 438, 0, 3, 96, 99,
                                                                       228, 234, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 448, 0, 3, 99,
                                                                       102, 234, 240, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 458, 0, 3, 102,
                                                                       105, 240, 246, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 468, 0, 3, 105,
                                                                       108, 246, 252, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 478, 0, 3, 114,
                                                                       120, 258, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 493, 0, 3, 120,
                                                                       126, 268, 278, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 508, 0, 3, 126,
                                                                       132, 278, 288, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 523, 0, 3, 132,
                                                                       138, 288, 298, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 538, 0, 3, 138,
                                                                       144, 298, 308, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 553, 0, 3, 144,
                                                                       150, 308, 318, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 568, 0, 3, 150,
                                                                       156, 318, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 583, 0, 3, 156,
                                                                       162, 328, 338, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 598, 0, 3, 162,
                                                                       168, 338, 348, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 613, 0, 3, 168,
                                                                       174, 348, 358, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 628, 0, 3, 186,
                                                                       192, 368, 378, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 643, 0, 3, 192,
                                                                       198, 378, 388, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 658, 0, 3, 198,
                                                                       204, 388, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 673, 0, 3, 204,
                                                                       210, 398, 408, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 688, 0, 3, 210,
                                                                       216, 408, 418, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 703, 0, 3, 216,
                                                                       222, 418, 428, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 718, 0, 3, 222,
                                                                       228, 428, 438, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 733, 0, 3, 228,
                                                                       234, 438, 448, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 748, 0, 3, 234,
                                                                       240, 448, 458, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 763, 0, 3, 240,
                                                                       246, 458, 468, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 778, 0, 3, 258,
                                                                       268, 478, 493, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 799, 0, 3, 268,
                                                                       278, 493, 508, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 820, 0, 3, 278,
                                                                       288, 508, 523, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 841, 0, 3, 288,
                                                                       298, 523, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 862, 0, 3, 298,
                                                                       308, 538, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 883, 0, 3, 308,
                                                                       318, 553, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 904, 0, 3, 318,
                                                                       328, 568, 583, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 925, 0, 3, 328,
                                                                       338, 583, 598, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 946, 0, 3, 338,
                                                                       348, 598, 613, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 967, 0, 3, 368,
                                                                       378, 628, 643, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 988, 0, 3, 378,
                                                                       388, 643, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1009, 0, 3, 388,
                                                                       398, 658, 673, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 398,
                                                                       408, 673, 688, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1051, 0, 3, 408,
                                                                       418, 688, 703, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 418,
                                                                       428, 703, 718, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1093, 0, 3, 428,
                                                                       438, 718, 733, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1114, 0, 3, 438,
                                                                       448, 733, 748, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 1135, 0, 3, 448,
                                                                       458, 748, 763, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 478,
                                                                       493, 778, 799, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1184, 0, 3, 493,
                                                                       508, 799, 820, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1212, 0, 3, 508,
                                                                       523, 820, 841, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1240, 0, 3, 523,
                                                                       538, 841, 862, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1268, 0, 3, 538,
                                                                       553, 862, 883, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1296, 0, 3, 553,
                                                                       568, 883, 904, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 568,
                                                                       583, 904, 925, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1352, 0, 3, 583,
                                                                       598, 925, 946, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1380, 0, 3, 628,
                                                                       643, 967, 988, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1408, 0, 3, 643,
                                                                       658, 988, 1009, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1436, 0, 3, 658,
                                                                       673, 1009, 1030, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1464, 0, 3, 673,
                                                                       688, 1030, 1051, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 688,
                                                                       703, 1051, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1520, 0, 3, 703,
                                                                       718, 1072, 1093, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1548, 0, 3, 718,
                                                                       733, 1093, 1114, ncols,
                                                                       gamma, p, q);

                    compute_prim_iss_three_center_electron_repulsion_0(buffer, 1576, 0, 3, 733,
                                                                       748, 1114, 1135, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1604, 0, 3, 778,
                                                                       799, 1156, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1640, 0, 3, 799,
                                                                       820, 1184, 1212, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1676, 0, 3, 820,
                                                                       841, 1212, 1240, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1712, 0, 3, 841,
                                                                       862, 1240, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1748, 0, 3, 862,
                                                                       883, 1268, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1784, 0, 3, 883,
                                                                       904, 1296, 1324, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1820, 0, 3, 904,
                                                                       925, 1324, 1352, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1856, 0, 3, 967,
                                                                       988, 1380, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1892, 0, 3, 988,
                                                                       1009, 1408, 1436, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1928, 0, 3, 1009,
                                                                       1030, 1436, 1464, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 1964, 0, 3, 1030,
                                                                       1051, 1464, 1492, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2000, 0, 3, 1051,
                                                                       1072, 1492, 1520, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2036, 0, 3, 1072,
                                                                       1093, 1520, 1548, ncols,
                                                                       gamma, p, q);

                    compute_prim_kss_three_center_electron_repulsion_0(buffer, 2072, 0, 3, 1093,
                                                                       1114, 1548, 1576, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2108, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2111, 3, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2114, 3, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2117, 3, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2120, 3, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2123, 3, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2126, 3, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2129, 3, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2132, 3, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2135, 3, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2138, 3, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2141, 3, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2144, 3, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2147, 3, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2150, 3, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2153, 3, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2156, 3, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2159, 3, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2162, 3, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2165, 3, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2168, 3, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2171, 3, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2174, 3, 34,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2177, 3, 35,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2180, 3, 9, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2189, 3, 10, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2198, 3, 11, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2207, 3, 12, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2216, 3, 13, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2225, 3, 14, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2234, 3, 15, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2243, 3, 16, 63,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2252, 3, 17, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2261, 3, 18, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2270, 3, 19, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2279, 3, 24, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2288, 3, 25, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2297, 3, 26, 87,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2306, 3, 27, 90,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2315, 3, 28, 93,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2324, 3, 29, 96,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2333, 3, 30, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2342, 3, 31, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2351, 3, 32, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2360, 3, 33, 108,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2369, 3, 34, 111,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2378, 3, 42, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2396, 3, 45, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2414, 3, 48, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2432, 3, 51, 144,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2450, 3, 54, 150,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2468, 3, 57, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2486, 3, 60, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2504, 3, 63, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2522, 3, 66, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2540, 3, 69, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2558, 3, 81, 198,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2576, 3, 84, 204,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2594, 3, 87, 210,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2612, 3, 90, 216,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2630, 3, 93, 222,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2648, 3, 96, 228,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2666, 3, 99, 234,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2684, 3, 102, 240,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2702, 3, 105, 246,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2720, 3, 108, 252,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2738, 3, 126, 278,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2768, 3, 132, 288,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2798, 3, 138, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2828, 3, 144, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2858, 3, 150, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2888, 3, 156, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2918, 3, 162, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2948, 3, 168, 348,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2978, 3, 174, 358,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3008, 3, 198, 388,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3038, 3, 204, 398,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3068, 3, 210, 408,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3098, 3, 216, 418,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3128, 3, 222, 428,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3158, 3, 228, 438,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3188, 3, 234, 448,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3218, 3, 240, 458,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3248, 3, 246, 468,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3278, 3, 278, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3323, 3, 288, 523,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3368, 3, 298, 538,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3413, 3, 308, 553,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3458, 3, 318, 568,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3503, 3, 328, 583,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3548, 3, 338, 598,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3593, 3, 348, 613,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3638, 3, 388, 658,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3683, 3, 398, 673,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3728, 3, 408, 688,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3773, 3, 418, 703,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3818, 3, 428, 718,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3863, 3, 438, 733,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3908, 3, 448, 748,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3953, 3, 458, 763,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 3998, 3, 508, 820,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4061, 3, 523, 841,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4124, 3, 538, 862,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4187, 3, 553, 883,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4250, 3, 568, 904,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4313, 3, 583, 925,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4376, 3, 598, 946,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4439, 3, 658,
                                                                       1009, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4502, 3, 673,
                                                                       1030, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4565, 3, 688,
                                                                       1051, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4628, 3, 703,
                                                                       1072, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4691, 3, 718,
                                                                       1093, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4754, 3, 733,
                                                                       1114, ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 4817, 3, 748,
                                                                       1135, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4880, 3, 820,
                                                                       1212, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 4964, 3, 841,
                                                                       1240, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5048, 3, 862,
                                                                       1268, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5132, 3, 883,
                                                                       1296, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5216, 3, 904,
                                                                       1324, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5300, 3, 925,
                                                                       1352, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5384, 3, 1009,
                                                                       1436, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5468, 3, 1030,
                                                                       1464, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5552, 3, 1051,
                                                                       1492, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5636, 3, 1072,
                                                                       1520, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5720, 3, 1093,
                                                                       1548, ncols, p, q);

                    compute_prim_isp_three_center_electron_repulsion_0(buffer, 5804, 3, 1114,
                                                                       1576, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5888, 3, 1212,
                                                                       1676, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 5996, 3, 1240,
                                                                       1712, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6104, 3, 1268,
                                                                       1748, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6212, 3, 1296,
                                                                       1784, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6320, 3, 1324,
                                                                       1820, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6428, 3, 1436,
                                                                       1928, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6536, 3, 1464,
                                                                       1964, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6644, 3, 1492,
                                                                       2000, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6752, 3, 1520,
                                                                       2036, ncols, p, q);

                    compute_prim_ksp_three_center_electron_repulsion_0(buffer, 6860, 3, 1548,
                                                                       2072, ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6968, 3, 7, 8,
                                                                       2108, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6974, 3, 8, 9,
                                                                       2111, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6980, 3, 9, 10,
                                                                       2114, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6986, 3, 10, 11,
                                                                       2117, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6992, 3, 11, 12,
                                                                       2120, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6998, 3, 12, 13,
                                                                       2123, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7004, 3, 13, 14,
                                                                       2126, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7010, 3, 14, 15,
                                                                       2129, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7016, 3, 15, 16,
                                                                       2132, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7022, 3, 16, 17,
                                                                       2135, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7028, 3, 17, 18,
                                                                       2138, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7034, 3, 18, 19,
                                                                       2141, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7040, 3, 22, 23,
                                                                       2144, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7046, 3, 23, 24,
                                                                       2147, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7052, 3, 24, 25,
                                                                       2150, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7058, 3, 25, 26,
                                                                       2153, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7064, 3, 26, 27,
                                                                       2156, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7070, 3, 27, 28,
                                                                       2159, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7076, 3, 28, 29,
                                                                       2162, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7082, 3, 29, 30,
                                                                       2165, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7088, 3, 30, 31,
                                                                       2168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7094, 3, 31, 32,
                                                                       2171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7100, 3, 32, 33,
                                                                       2174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7106, 3, 33, 34,
                                                                       2177, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7112, 0, 3, 6968,
                                                                       2108, 6974, 2180, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7130, 0, 3, 6974,
                                                                       2111, 6980, 2189, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7148, 0, 3, 6980,
                                                                       2114, 6986, 2198, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7166, 0, 3, 6986,
                                                                       2117, 6992, 2207, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7184, 0, 3, 6992,
                                                                       2120, 6998, 2216, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7202, 0, 3, 6998,
                                                                       2123, 7004, 2225, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7220, 0, 3, 7004,
                                                                       2126, 7010, 2234, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7238, 0, 3, 7010,
                                                                       2129, 7016, 2243, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7256, 0, 3, 7016,
                                                                       2132, 7022, 2252, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7274, 0, 3, 7022,
                                                                       2135, 7028, 2261, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7292, 0, 3, 7028,
                                                                       2138, 7034, 2270, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7310, 0, 3, 7040,
                                                                       2144, 7046, 2279, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7328, 0, 3, 7046,
                                                                       2147, 7052, 2288, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7346, 0, 3, 7052,
                                                                       2150, 7058, 2297, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7364, 0, 3, 7058,
                                                                       2153, 7064, 2306, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7382, 0, 3, 7064,
                                                                       2156, 7070, 2315, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7400, 0, 3, 7070,
                                                                       2159, 7076, 2324, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7418, 0, 3, 7076,
                                                                       2162, 7082, 2333, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7436, 0, 3, 7082,
                                                                       2165, 7088, 2342, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7454, 0, 3, 7088,
                                                                       2168, 7094, 2351, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7472, 0, 3, 7094,
                                                                       2171, 7100, 2360, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7490, 0, 3, 7100,
                                                                       2174, 7106, 2369, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7508, 0, 3, 7112,
                                                                       2180, 7130, 114, 120,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7544, 0, 3, 7130,
                                                                       2189, 7148, 120, 126,
                                                                       2396, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7580, 0, 3, 7148,
                                                                       2198, 7166, 126, 132,
                                                                       2414, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7616, 0, 3, 7166,
                                                                       2207, 7184, 132, 138,
                                                                       2432, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7652, 0, 3, 7184,
                                                                       2216, 7202, 138, 144,
                                                                       2450, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7688, 0, 3, 7202,
                                                                       2225, 7220, 144, 150,
                                                                       2468, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7724, 0, 3, 7220,
                                                                       2234, 7238, 150, 156,
                                                                       2486, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7760, 0, 3, 7238,
                                                                       2243, 7256, 156, 162,
                                                                       2504, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7796, 0, 3, 7256,
                                                                       2252, 7274, 162, 168,
                                                                       2522, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7832, 0, 3, 7274,
                                                                       2261, 7292, 168, 174,
                                                                       2540, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7868, 0, 3, 7310,
                                                                       2279, 7328, 186, 192,
                                                                       2558, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7904, 0, 3, 7328,
                                                                       2288, 7346, 192, 198,
                                                                       2576, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7940, 0, 3, 7346,
                                                                       2297, 7364, 198, 204,
                                                                       2594, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 7976, 0, 3, 7364,
                                                                       2306, 7382, 204, 210,
                                                                       2612, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8012, 0, 3, 7382,
                                                                       2315, 7400, 210, 216,
                                                                       2630, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8048, 0, 3, 7400,
                                                                       2324, 7418, 216, 222,
                                                                       2648, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8084, 0, 3, 7418,
                                                                       2333, 7436, 222, 228,
                                                                       2666, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8120, 0, 3, 7436,
                                                                       2342, 7454, 228, 234,
                                                                       2684, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8156, 0, 3, 7454,
                                                                       2351, 7472, 234, 240,
                                                                       2702, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8192, 0, 3, 7472,
                                                                       2360, 7490, 240, 246,
                                                                       2720, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8228, 0, 3, 7508,
                                                                       2378, 7544, 258, 268,
                                                                       2738, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8288, 0, 3, 7544,
                                                                       2396, 7580, 268, 278,
                                                                       2768, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8348, 0, 3, 7580,
                                                                       2414, 7616, 278, 288,
                                                                       2798, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8408, 0, 3, 7616,
                                                                       2432, 7652, 288, 298,
                                                                       2828, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8468, 0, 3, 7652,
                                                                       2450, 7688, 298, 308,
                                                                       2858, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8528, 0, 3, 7688,
                                                                       2468, 7724, 308, 318,
                                                                       2888, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8588, 0, 3, 7724,
                                                                       2486, 7760, 318, 328,
                                                                       2918, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8648, 0, 3, 7760,
                                                                       2504, 7796, 328, 338,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8708, 0, 3, 7796,
                                                                       2522, 7832, 338, 348,
                                                                       2978, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8768, 0, 3, 7868,
                                                                       2558, 7904, 368, 378,
                                                                       3008, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8828, 0, 3, 7904,
                                                                       2576, 7940, 378, 388,
                                                                       3038, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8888, 0, 3, 7940,
                                                                       2594, 7976, 388, 398,
                                                                       3068, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8948, 0, 3, 7976,
                                                                       2612, 8012, 398, 408,
                                                                       3098, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9008, 0, 3, 8012,
                                                                       2630, 8048, 408, 418,
                                                                       3128, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9068, 0, 3, 8048,
                                                                       2648, 8084, 418, 428,
                                                                       3158, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9128, 0, 3, 8084,
                                                                       2666, 8120, 428, 438,
                                                                       3188, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9188, 0, 3, 8120,
                                                                       2684, 8156, 438, 448,
                                                                       3218, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9248, 0, 3, 8156,
                                                                       2702, 8192, 448, 458,
                                                                       3248, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9308, 0, 3, 8228,
                                                                       2738, 8288, 478, 493,
                                                                       3278, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9398, 0, 3, 8288,
                                                                       2768, 8348, 493, 508,
                                                                       3323, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9488, 0, 3, 8348,
                                                                       2798, 8408, 508, 523,
                                                                       3368, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9578, 0, 3, 8408,
                                                                       2828, 8468, 523, 538,
                                                                       3413, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9668, 0, 3, 8468,
                                                                       2858, 8528, 538, 553,
                                                                       3458, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9758, 0, 3, 8528,
                                                                       2888, 8588, 553, 568,
                                                                       3503, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9848, 0, 3, 8588,
                                                                       2918, 8648, 568, 583,
                                                                       3548, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 9938, 0, 3, 8648,
                                                                       2948, 8708, 583, 598,
                                                                       3593, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10028, 0, 3, 8768,
                                                                       3008, 8828, 628, 643,
                                                                       3638, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10118, 0, 3, 8828,
                                                                       3038, 8888, 643, 658,
                                                                       3683, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10208, 0, 3, 8888,
                                                                       3068, 8948, 658, 673,
                                                                       3728, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10298, 0, 3, 8948,
                                                                       3098, 9008, 673, 688,
                                                                       3773, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10388, 0, 3, 9008,
                                                                       3128, 9068, 688, 703,
                                                                       3818, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10478, 0, 3, 9068,
                                                                       3158, 9128, 703, 718,
                                                                       3863, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10568, 0, 3, 9128,
                                                                       3188, 9188, 718, 733,
                                                                       3908, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10658, 0, 3, 9188,
                                                                       3218, 9248, 733, 748,
                                                                       3953, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10748, 0, 3, 9308,
                                                                       3278, 9398, 778, 799,
                                                                       3998, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 10874, 0, 3, 9398,
                                                                       3323, 9488, 799, 820,
                                                                       4061, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11000, 0, 3, 9488,
                                                                       3368, 9578, 820, 841,
                                                                       4124, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11126, 0, 3, 9578,
                                                                       3413, 9668, 841, 862,
                                                                       4187, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11252, 0, 3, 9668,
                                                                       3458, 9758, 862, 883,
                                                                       4250, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11378, 0, 3, 9758,
                                                                       3503, 9848, 883, 904,
                                                                       4313, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11504, 0, 3, 9848,
                                                                       3548, 9938, 904, 925,
                                                                       4376, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11630, 0, 3,
                                                                       10028, 3638, 10118, 967,
                                                                       988, 4439, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11756, 0, 3,
                                                                       10118, 3683, 10208, 988,
                                                                       1009, 4502, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 11882, 0, 3,
                                                                       10208, 3728, 10298, 1009,
                                                                       1030, 4565, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12008, 0, 3,
                                                                       10298, 3773, 10388, 1030,
                                                                       1051, 4628, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12134, 0, 3,
                                                                       10388, 3818, 10478, 1051,
                                                                       1072, 4691, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12260, 0, 3,
                                                                       10478, 3863, 10568, 1072,
                                                                       1093, 4754, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 12386, 0, 3,
                                                                       10568, 3908, 10658, 1093,
                                                                       1114, 4817, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12512, 0, 3,
                                                                       10748, 3998, 10874, 1156,
                                                                       1184, 4880, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12680, 0, 3,
                                                                       10874, 4061, 11000, 1184,
                                                                       1212, 4964, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 12848, 0, 3,
                                                                       11000, 4124, 11126, 1212,
                                                                       1240, 5048, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13016, 0, 3,
                                                                       11126, 4187, 11252, 1240,
                                                                       1268, 5132, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13184, 0, 3,
                                                                       11252, 4250, 11378, 1268,
                                                                       1296, 5216, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13352, 0, 3,
                                                                       11378, 4313, 11504, 1296,
                                                                       1324, 5300, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13520, 0, 3,
                                                                       11630, 4439, 11756, 1380,
                                                                       1408, 5384, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13688, 0, 3,
                                                                       11756, 4502, 11882, 1408,
                                                                       1436, 5468, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 13856, 0, 3,
                                                                       11882, 4565, 12008, 1436,
                                                                       1464, 5552, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14024, 0, 3,
                                                                       12008, 4628, 12134, 1464,
                                                                       1492, 5636, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14192, 0, 3,
                                                                       12134, 4691, 12260, 1492,
                                                                       1520, 5720, ncols, gamma,
                                                                       p, q);

                    compute_prim_isd_three_center_electron_repulsion_0(buffer, 14360, 0, 3,
                                                                       12260, 4754, 12386, 1520,
                                                                       1548, 5804, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14528, 0, 3,
                                                                       12512, 4880, 12680, 1604,
                                                                       1640, 5888, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14744, 0, 3,
                                                                       12680, 4964, 12848, 1640,
                                                                       1676, 5996, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 14960, 0, 3,
                                                                       12848, 5048, 13016, 1676,
                                                                       1712, 6104, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15176, 0, 3,
                                                                       13016, 5132, 13184, 1712,
                                                                       1748, 6212, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15392, 0, 3,
                                                                       13184, 5216, 13352, 1748,
                                                                       1784, 6320, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15608, 0, 3,
                                                                       13520, 5384, 13688, 1856,
                                                                       1892, 6428, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 15824, 0, 3,
                                                                       13688, 5468, 13856, 1892,
                                                                       1928, 6536, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16040, 0, 3,
                                                                       13856, 5552, 14024, 1928,
                                                                       1964, 6644, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16256, 0, 3,
                                                                       14024, 5636, 14192, 1964,
                                                                       2000, 6752, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksd_three_center_electron_repulsion_0(buffer, 16472, 0, 3,
                                                                       14192, 5720, 14360, 2000,
                                                                       2036, 6860, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16688, 3, 2108,
                                                                       2111, 6980, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16698, 3, 2111,
                                                                       2114, 6986, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16708, 3, 2114,
                                                                       2117, 6992, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16718, 3, 2117,
                                                                       2120, 6998, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16728, 3, 2120,
                                                                       2123, 7004, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16738, 3, 2123,
                                                                       2126, 7010, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16748, 3, 2126,
                                                                       2129, 7016, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16758, 3, 2129,
                                                                       2132, 7022, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16768, 3, 2132,
                                                                       2135, 7028, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16778, 3, 2135,
                                                                       2138, 7034, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16788, 3, 2144,
                                                                       2147, 7052, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16798, 3, 2147,
                                                                       2150, 7058, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16808, 3, 2150,
                                                                       2153, 7064, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16818, 3, 2153,
                                                                       2156, 7070, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16828, 3, 2156,
                                                                       2159, 7076, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16838, 3, 2159,
                                                                       2162, 7082, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16848, 3, 2162,
                                                                       2165, 7088, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16858, 3, 2165,
                                                                       2168, 7094, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16868, 3, 2168,
                                                                       2171, 7100, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 16878, 3, 2171,
                                                                       2174, 7106, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16888, 0, 3,
                                                                       16688, 6980, 16698, 2180,
                                                                       2189, 7148, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16918, 0, 3,
                                                                       16698, 6986, 16708, 2189,
                                                                       2198, 7166, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16948, 0, 3,
                                                                       16708, 6992, 16718, 2198,
                                                                       2207, 7184, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 16978, 0, 3,
                                                                       16718, 6998, 16728, 2207,
                                                                       2216, 7202, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17008, 0, 3,
                                                                       16728, 7004, 16738, 2216,
                                                                       2225, 7220, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17038, 0, 3,
                                                                       16738, 7010, 16748, 2225,
                                                                       2234, 7238, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17068, 0, 3,
                                                                       16748, 7016, 16758, 2234,
                                                                       2243, 7256, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17098, 0, 3,
                                                                       16758, 7022, 16768, 2243,
                                                                       2252, 7274, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17128, 0, 3,
                                                                       16768, 7028, 16778, 2252,
                                                                       2261, 7292, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17158, 0, 3,
                                                                       16788, 7052, 16798, 2279,
                                                                       2288, 7346, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17188, 0, 3,
                                                                       16798, 7058, 16808, 2288,
                                                                       2297, 7364, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17218, 0, 3,
                                                                       16808, 7064, 16818, 2297,
                                                                       2306, 7382, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17248, 0, 3,
                                                                       16818, 7070, 16828, 2306,
                                                                       2315, 7400, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17278, 0, 3,
                                                                       16828, 7076, 16838, 2315,
                                                                       2324, 7418, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17308, 0, 3,
                                                                       16838, 7082, 16848, 2324,
                                                                       2333, 7436, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17338, 0, 3,
                                                                       16848, 7088, 16858, 2333,
                                                                       2342, 7454, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17368, 0, 3,
                                                                       16858, 7094, 16868, 2342,
                                                                       2351, 7472, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17398, 0, 3,
                                                                       16868, 7100, 16878, 2351,
                                                                       2360, 7490, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17428, 0, 3,
                                                                       16888, 7148, 16918, 2378,
                                                                       2396, 7580, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17488, 0, 3,
                                                                       16918, 7166, 16948, 2396,
                                                                       2414, 7616, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17548, 0, 3,
                                                                       16948, 7184, 16978, 2414,
                                                                       2432, 7652, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17608, 0, 3,
                                                                       16978, 7202, 17008, 2432,
                                                                       2450, 7688, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17668, 0, 3,
                                                                       17008, 7220, 17038, 2450,
                                                                       2468, 7724, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17728, 0, 3,
                                                                       17038, 7238, 17068, 2468,
                                                                       2486, 7760, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17788, 0, 3,
                                                                       17068, 7256, 17098, 2486,
                                                                       2504, 7796, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17848, 0, 3,
                                                                       17098, 7274, 17128, 2504,
                                                                       2522, 7832, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17908, 0, 3,
                                                                       17158, 7346, 17188, 2558,
                                                                       2576, 7940, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 17968, 0, 3,
                                                                       17188, 7364, 17218, 2576,
                                                                       2594, 7976, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18028, 0, 3,
                                                                       17218, 7382, 17248, 2594,
                                                                       2612, 8012, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18088, 0, 3,
                                                                       17248, 7400, 17278, 2612,
                                                                       2630, 8048, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18148, 0, 3,
                                                                       17278, 7418, 17308, 2630,
                                                                       2648, 8084, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18208, 0, 3,
                                                                       17308, 7436, 17338, 2648,
                                                                       2666, 8120, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18268, 0, 3,
                                                                       17338, 7454, 17368, 2666,
                                                                       2684, 8156, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 18328, 0, 3,
                                                                       17368, 7472, 17398, 2684,
                                                                       2702, 8192, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18388, 0, 3,
                                                                       17428, 7580, 17488, 2738,
                                                                       2768, 8348, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18488, 0, 3,
                                                                       17488, 7616, 17548, 2768,
                                                                       2798, 8408, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18588, 0, 3,
                                                                       17548, 7652, 17608, 2798,
                                                                       2828, 8468, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18688, 0, 3,
                                                                       17608, 7688, 17668, 2828,
                                                                       2858, 8528, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18788, 0, 3,
                                                                       17668, 7724, 17728, 2858,
                                                                       2888, 8588, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18888, 0, 3,
                                                                       17728, 7760, 17788, 2888,
                                                                       2918, 8648, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18988, 0, 3,
                                                                       17788, 7796, 17848, 2918,
                                                                       2948, 8708, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19088, 0, 3,
                                                                       17908, 7940, 17968, 3008,
                                                                       3038, 8888, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19188, 0, 3,
                                                                       17968, 7976, 18028, 3038,
                                                                       3068, 8948, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19288, 0, 3,
                                                                       18028, 8012, 18088, 3068,
                                                                       3098, 9008, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19388, 0, 3,
                                                                       18088, 8048, 18148, 3098,
                                                                       3128, 9068, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19488, 0, 3,
                                                                       18148, 8084, 18208, 3128,
                                                                       3158, 9128, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19588, 0, 3,
                                                                       18208, 8120, 18268, 3158,
                                                                       3188, 9188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19688, 0, 3,
                                                                       18268, 8156, 18328, 3188,
                                                                       3218, 9248, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 19788, 0, 3,
                                                                       18388, 8348, 18488, 3278,
                                                                       3323, 9488, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 19938, 0, 3,
                                                                       18488, 8408, 18588, 3323,
                                                                       3368, 9578, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20088, 0, 3,
                                                                       18588, 8468, 18688, 3368,
                                                                       3413, 9668, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20238, 0, 3,
                                                                       18688, 8528, 18788, 3413,
                                                                       3458, 9758, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20388, 0, 3,
                                                                       18788, 8588, 18888, 3458,
                                                                       3503, 9848, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20538, 0, 3,
                                                                       18888, 8648, 18988, 3503,
                                                                       3548, 9938, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20688, 0, 3,
                                                                       19088, 8888, 19188, 3638,
                                                                       3683, 10208, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20838, 0, 3,
                                                                       19188, 8948, 19288, 3683,
                                                                       3728, 10298, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 20988, 0, 3,
                                                                       19288, 9008, 19388, 3728,
                                                                       3773, 10388, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21138, 0, 3,
                                                                       19388, 9068, 19488, 3773,
                                                                       3818, 10478, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21288, 0, 3,
                                                                       19488, 9128, 19588, 3818,
                                                                       3863, 10568, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 21438, 0, 3,
                                                                       19588, 9188, 19688, 3863,
                                                                       3908, 10658, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 21588, 0, 3,
                                                                       19788, 9488, 19938, 3998,
                                                                       4061, 11000, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 21798, 0, 3,
                                                                       19938, 9578, 20088, 4061,
                                                                       4124, 11126, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22008, 0, 3,
                                                                       20088, 9668, 20238, 4124,
                                                                       4187, 11252, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22218, 0, 3,
                                                                       20238, 9758, 20388, 4187,
                                                                       4250, 11378, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22428, 0, 3,
                                                                       20388, 9848, 20538, 4250,
                                                                       4313, 11504, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22638, 0, 3,
                                                                       20688, 10208, 20838, 4439,
                                                                       4502, 11882, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 22848, 0, 3,
                                                                       20838, 10298, 20988, 4502,
                                                                       4565, 12008, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 23058, 0, 3,
                                                                       20988, 10388, 21138, 4565,
                                                                       4628, 12134, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 23268, 0, 3,
                                                                       21138, 10478, 21288, 4628,
                                                                       4691, 12260, ncols, gamma,
                                                                       p, q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 23478, 0, 3,
                                                                       21288, 10568, 21438, 4691,
                                                                       4754, 12386, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 23688, 0, 3,
                                                                       21588, 11000, 21798, 4880,
                                                                       4964, 12848, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 23968, 0, 3,
                                                                       21798, 11126, 22008, 4964,
                                                                       5048, 13016, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24248, 0, 3,
                                                                       22008, 11252, 22218, 5048,
                                                                       5132, 13184, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24528, 0, 3,
                                                                       22218, 11378, 22428, 5132,
                                                                       5216, 13352, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 24808, 0, 3,
                                                                       22638, 11882, 22848, 5384,
                                                                       5468, 13856, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 25088, 0, 3,
                                                                       22848, 12008, 23058, 5468,
                                                                       5552, 14024, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 25368, 0, 3,
                                                                       23058, 12134, 23268, 5552,
                                                                       5636, 14192, ncols, gamma,
                                                                       p, q);

                    compute_prim_isf_three_center_electron_repulsion_0(buffer, 25648, 0, 3,
                                                                       23268, 12260, 23478, 5636,
                                                                       5720, 14360, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 25928, 0, 3,
                                                                       23688, 12848, 23968, 5888,
                                                                       5996, 14960, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 26288, 0, 3,
                                                                       23968, 13016, 24248, 5996,
                                                                       6104, 15176, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 26648, 0, 3,
                                                                       24248, 13184, 24528, 6104,
                                                                       6212, 15392, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 27008, 0, 3,
                                                                       24808, 13856, 25088, 6428,
                                                                       6536, 16040, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 27368, 0, 3,
                                                                       25088, 14024, 25368, 6536,
                                                                       6644, 16256, ncols, gamma,
                                                                       p, q);

                    compute_prim_ksf_three_center_electron_repulsion_0(buffer, 27728, 0, 3,
                                                                       25368, 14192, 25648, 6644,
                                                                       6752, 16472, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28088, 3, 6968,
                                                                       6974, 16688, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28103, 3, 6974,
                                                                       6980, 16698, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28118, 3, 6980,
                                                                       6986, 16708, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28133, 3, 6986,
                                                                       6992, 16718, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28148, 3, 6992,
                                                                       6998, 16728, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28163, 3, 6998,
                                                                       7004, 16738, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28178, 3, 7004,
                                                                       7010, 16748, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28193, 3, 7010,
                                                                       7016, 16758, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28208, 3, 7016,
                                                                       7022, 16768, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28223, 3, 7022,
                                                                       7028, 16778, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28238, 3, 7040,
                                                                       7046, 16788, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28253, 3, 7046,
                                                                       7052, 16798, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28268, 3, 7052,
                                                                       7058, 16808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28283, 3, 7058,
                                                                       7064, 16818, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28298, 3, 7064,
                                                                       7070, 16828, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28313, 3, 7070,
                                                                       7076, 16838, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28328, 3, 7076,
                                                                       7082, 16848, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28343, 3, 7082,
                                                                       7088, 16858, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28358, 3, 7088,
                                                                       7094, 16868, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 28373, 3, 7094,
                                                                       7100, 16878, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28388, 0, 3,
                                                                       28088, 16688, 28103, 7112,
                                                                       7130, 16888, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28433, 0, 3,
                                                                       28103, 16698, 28118, 7130,
                                                                       7148, 16918, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28478, 0, 3,
                                                                       28118, 16708, 28133, 7148,
                                                                       7166, 16948, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28523, 0, 3,
                                                                       28133, 16718, 28148, 7166,
                                                                       7184, 16978, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28568, 0, 3,
                                                                       28148, 16728, 28163, 7184,
                                                                       7202, 17008, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28613, 0, 3,
                                                                       28163, 16738, 28178, 7202,
                                                                       7220, 17038, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28658, 0, 3,
                                                                       28178, 16748, 28193, 7220,
                                                                       7238, 17068, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28703, 0, 3,
                                                                       28193, 16758, 28208, 7238,
                                                                       7256, 17098, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28748, 0, 3,
                                                                       28208, 16768, 28223, 7256,
                                                                       7274, 17128, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28793, 0, 3,
                                                                       28238, 16788, 28253, 7310,
                                                                       7328, 17158, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28838, 0, 3,
                                                                       28253, 16798, 28268, 7328,
                                                                       7346, 17188, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28883, 0, 3,
                                                                       28268, 16808, 28283, 7346,
                                                                       7364, 17218, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28928, 0, 3,
                                                                       28283, 16818, 28298, 7364,
                                                                       7382, 17248, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 28973, 0, 3,
                                                                       28298, 16828, 28313, 7382,
                                                                       7400, 17278, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29018, 0, 3,
                                                                       28313, 16838, 28328, 7400,
                                                                       7418, 17308, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29063, 0, 3,
                                                                       28328, 16848, 28343, 7418,
                                                                       7436, 17338, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29108, 0, 3,
                                                                       28343, 16858, 28358, 7436,
                                                                       7454, 17368, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29153, 0, 3,
                                                                       28358, 16868, 28373, 7454,
                                                                       7472, 17398, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 29198, 0, 3,
                                                                       28388, 16888, 28433, 7508,
                                                                       7544, 17428, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 29288, 0, 3,
                                                                       28433, 16918, 28478, 7544,
                                                                       7580, 17488, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 29378, 0, 3,
                                                                       28478, 16948, 28523, 7580,
                                                                       7616, 17548, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 29468, 0, 3,
                                                                       28523, 16978, 28568, 7616,
                                                                       7652, 17608, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 29558, 0, 3,
                                                                       28568, 17008, 28613, 7652,
                                                                       7688, 17668, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 29648, 0, 3,
                                                                       28613, 17038, 28658, 7688,
                                                                       7724, 17728, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 29738, 0, 3,
                                                                       28658, 17068, 28703, 7724,
                                                                       7760, 17788, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 29828, 0, 3,
                                                                       28703, 17098, 28748, 7760,
                                                                       7796, 17848, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 29918, 0, 3,
                                                                       28793, 17158, 28838, 7868,
                                                                       7904, 17908, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 30008, 0, 3,
                                                                       28838, 17188, 28883, 7904,
                                                                       7940, 17968, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 30098, 0, 3,
                                                                       28883, 17218, 28928, 7940,
                                                                       7976, 18028, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 30188, 0, 3,
                                                                       28928, 17248, 28973, 7976,
                                                                       8012, 18088, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 30278, 0, 3,
                                                                       28973, 17278, 29018, 8012,
                                                                       8048, 18148, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 30368, 0, 3,
                                                                       29018, 17308, 29063, 8048,
                                                                       8084, 18208, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 30458, 0, 3,
                                                                       29063, 17338, 29108, 8084,
                                                                       8120, 18268, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 30548, 0, 3,
                                                                       29108, 17368, 29153, 8120,
                                                                       8156, 18328, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 30638, 0, 3,
                                                                       29198, 17428, 29288, 8228,
                                                                       8288, 18388, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 30788, 0, 3,
                                                                       29288, 17488, 29378, 8288,
                                                                       8348, 18488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 30938, 0, 3,
                                                                       29378, 17548, 29468, 8348,
                                                                       8408, 18588, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 31088, 0, 3,
                                                                       29468, 17608, 29558, 8408,
                                                                       8468, 18688, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 31238, 0, 3,
                                                                       29558, 17668, 29648, 8468,
                                                                       8528, 18788, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 31388, 0, 3,
                                                                       29648, 17728, 29738, 8528,
                                                                       8588, 18888, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 31538, 0, 3,
                                                                       29738, 17788, 29828, 8588,
                                                                       8648, 18988, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 31688, 0, 3,
                                                                       29918, 17908, 30008, 8768,
                                                                       8828, 19088, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 31838, 0, 3,
                                                                       30008, 17968, 30098, 8828,
                                                                       8888, 19188, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 31988, 0, 3,
                                                                       30098, 18028, 30188, 8888,
                                                                       8948, 19288, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 32138, 0, 3,
                                                                       30188, 18088, 30278, 8948,
                                                                       9008, 19388, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 32288, 0, 3,
                                                                       30278, 18148, 30368, 9008,
                                                                       9068, 19488, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 32438, 0, 3,
                                                                       30368, 18208, 30458, 9068,
                                                                       9128, 19588, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 32588, 0, 3,
                                                                       30458, 18268, 30548, 9128,
                                                                       9188, 19688, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 32738, 0, 3,
                                                                       30638, 18388, 30788, 9308,
                                                                       9398, 19788, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 32963, 0, 3,
                                                                       30788, 18488, 30938, 9398,
                                                                       9488, 19938, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 33188, 0, 3,
                                                                       30938, 18588, 31088, 9488,
                                                                       9578, 20088, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 33413, 0, 3,
                                                                       31088, 18688, 31238, 9578,
                                                                       9668, 20238, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 33638, 0, 3,
                                                                       31238, 18788, 31388, 9668,
                                                                       9758, 20388, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 33863, 0, 3,
                                                                       31388, 18888, 31538, 9758,
                                                                       9848, 20538, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 34088, 0, 3,
                                                                       31688, 19088, 31838,
                                                                       10028, 10118, 20688,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 34313, 0, 3,
                                                                       31838, 19188, 31988,
                                                                       10118, 10208, 20838,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 34538, 0, 3,
                                                                       31988, 19288, 32138,
                                                                       10208, 10298, 20988,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 34763, 0, 3,
                                                                       32138, 19388, 32288,
                                                                       10298, 10388, 21138,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 34988, 0, 3,
                                                                       32288, 19488, 32438,
                                                                       10388, 10478, 21288,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 35213, 0, 3,
                                                                       32438, 19588, 32588,
                                                                       10478, 10568, 21438,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 35438, 0, 3,
                                                                       32738, 19788, 32963,
                                                                       10748, 10874, 21588,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 35753, 0, 3,
                                                                       32963, 19938, 33188,
                                                                       10874, 11000, 21798,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 36068, 0, 3,
                                                                       33188, 20088, 33413,
                                                                       11000, 11126, 22008,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 36383, 0, 3,
                                                                       33413, 20238, 33638,
                                                                       11126, 11252, 22218,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 36698, 0, 3,
                                                                       33638, 20388, 33863,
                                                                       11252, 11378, 22428,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 37013, 0, 3,
                                                                       34088, 20688, 34313,
                                                                       11630, 11756, 22638,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 37328, 0, 3,
                                                                       34313, 20838, 34538,
                                                                       11756, 11882, 22848,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 37643, 0, 3,
                                                                       34538, 20988, 34763,
                                                                       11882, 12008, 23058,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 37958, 0, 3,
                                                                       34763, 21138, 34988,
                                                                       12008, 12134, 23268,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsg_three_center_electron_repulsion_0(buffer, 38273, 0, 3,
                                                                       34988, 21288, 35213,
                                                                       12134, 12260, 23478,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 38588, 0, 3,
                                                                       35438, 21588, 35753,
                                                                       12512, 12680, 23688,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 39008, 0, 3,
                                                                       35753, 21798, 36068,
                                                                       12680, 12848, 23968,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 39428, 0, 3,
                                                                       36068, 22008, 36383,
                                                                       12848, 13016, 24248,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 39848, 0, 3,
                                                                       36383, 22218, 36698,
                                                                       13016, 13184, 24528,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 40268, 0, 3,
                                                                       37013, 22638, 37328,
                                                                       13520, 13688, 24808,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 40688, 0, 3,
                                                                       37328, 22848, 37643,
                                                                       13688, 13856, 25088,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 41108, 0, 3,
                                                                       37643, 23058, 37958,
                                                                       13856, 14024, 25368,
                                                                       ncols, gamma, p, q);

                    compute_prim_isg_three_center_electron_repulsion_0(buffer, 41528, 0, 3,
                                                                       37958, 23268, 38273,
                                                                       14024, 14192, 25648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 41948, 0, 3,
                                                                       38588, 23688, 39008,
                                                                       14528, 14744, 25928,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 42488, 0, 3,
                                                                       39008, 23968, 39428,
                                                                       14744, 14960, 26288,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 43028, 0, 3,
                                                                       39428, 24248, 39848,
                                                                       14960, 15176, 26648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 43568, 0, 3,
                                                                       40268, 24808, 40688,
                                                                       15608, 15824, 27008,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 44108, 0, 3,
                                                                       40688, 25088, 41108,
                                                                       15824, 16040, 27368,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksg_three_center_electron_repulsion_0(buffer, 44648, 0, 3,
                                                                       41108, 25368, 41528,
                                                                       16040, 16256, 27728,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45188, 3, 16688,
                                                                       16698, 28118, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45209, 3, 16698,
                                                                       16708, 28133, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45230, 3, 16708,
                                                                       16718, 28148, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45251, 3, 16718,
                                                                       16728, 28163, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45272, 3, 16728,
                                                                       16738, 28178, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45293, 3, 16738,
                                                                       16748, 28193, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45314, 3, 16748,
                                                                       16758, 28208, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45335, 3, 16758,
                                                                       16768, 28223, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45356, 3, 16788,
                                                                       16798, 28268, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45377, 3, 16798,
                                                                       16808, 28283, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45398, 3, 16808,
                                                                       16818, 28298, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45419, 3, 16818,
                                                                       16828, 28313, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45440, 3, 16828,
                                                                       16838, 28328, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45461, 3, 16838,
                                                                       16848, 28343, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45482, 3, 16848,
                                                                       16858, 28358, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 45503, 3, 16858,
                                                                       16868, 28373, ncols,
                                                                       gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 45524, 0, 3,
                                                                       45188, 28118, 45209,
                                                                       16888, 16918, 28478,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 45587, 0, 3,
                                                                       45209, 28133, 45230,
                                                                       16918, 16948, 28523,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 45650, 0, 3,
                                                                       45230, 28148, 45251,
                                                                       16948, 16978, 28568,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 45713, 0, 3,
                                                                       45251, 28163, 45272,
                                                                       16978, 17008, 28613,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 45776, 0, 3,
                                                                       45272, 28178, 45293,
                                                                       17008, 17038, 28658,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 45839, 0, 3,
                                                                       45293, 28193, 45314,
                                                                       17038, 17068, 28703,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 45902, 0, 3,
                                                                       45314, 28208, 45335,
                                                                       17068, 17098, 28748,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 45965, 0, 3,
                                                                       45356, 28268, 45377,
                                                                       17158, 17188, 28883,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 46028, 0, 3,
                                                                       45377, 28283, 45398,
                                                                       17188, 17218, 28928,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 46091, 0, 3,
                                                                       45398, 28298, 45419,
                                                                       17218, 17248, 28973,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 46154, 0, 3,
                                                                       45419, 28313, 45440,
                                                                       17248, 17278, 29018,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 46217, 0, 3,
                                                                       45440, 28328, 45461,
                                                                       17278, 17308, 29063,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 46280, 0, 3,
                                                                       45461, 28343, 45482,
                                                                       17308, 17338, 29108,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 46343, 0, 3,
                                                                       45482, 28358, 45503,
                                                                       17338, 17368, 29153,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 46406, 0, 3,
                                                                       45524, 28478, 45587,
                                                                       17428, 17488, 29378,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 46532, 0, 3,
                                                                       45587, 28523, 45650,
                                                                       17488, 17548, 29468,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 46658, 0, 3,
                                                                       45650, 28568, 45713,
                                                                       17548, 17608, 29558,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 46784, 0, 3,
                                                                       45713, 28613, 45776,
                                                                       17608, 17668, 29648,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 46910, 0, 3,
                                                                       45776, 28658, 45839,
                                                                       17668, 17728, 29738,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 47036, 0, 3,
                                                                       45839, 28703, 45902,
                                                                       17728, 17788, 29828,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 47162, 0, 3,
                                                                       45965, 28883, 46028,
                                                                       17908, 17968, 30098,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 47288, 0, 3,
                                                                       46028, 28928, 46091,
                                                                       17968, 18028, 30188,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 47414, 0, 3,
                                                                       46091, 28973, 46154,
                                                                       18028, 18088, 30278,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 47540, 0, 3,
                                                                       46154, 29018, 46217,
                                                                       18088, 18148, 30368,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 47666, 0, 3,
                                                                       46217, 29063, 46280,
                                                                       18148, 18208, 30458,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 47792, 0, 3,
                                                                       46280, 29108, 46343,
                                                                       18208, 18268, 30548,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 47918, 0, 3,
                                                                       46406, 29378, 46532,
                                                                       18388, 18488, 30938,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 48128, 0, 3,
                                                                       46532, 29468, 46658,
                                                                       18488, 18588, 31088,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 48338, 0, 3,
                                                                       46658, 29558, 46784,
                                                                       18588, 18688, 31238,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 48548, 0, 3,
                                                                       46784, 29648, 46910,
                                                                       18688, 18788, 31388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 48758, 0, 3,
                                                                       46910, 29738, 47036,
                                                                       18788, 18888, 31538,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 48968, 0, 3,
                                                                       47162, 30098, 47288,
                                                                       19088, 19188, 31988,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 49178, 0, 3,
                                                                       47288, 30188, 47414,
                                                                       19188, 19288, 32138,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 49388, 0, 3,
                                                                       47414, 30278, 47540,
                                                                       19288, 19388, 32288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 49598, 0, 3,
                                                                       47540, 30368, 47666,
                                                                       19388, 19488, 32438,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 49808, 0, 3,
                                                                       47666, 30458, 47792,
                                                                       19488, 19588, 32588,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 50018, 0, 3,
                                                                       47918, 30938, 48128,
                                                                       19788, 19938, 33188,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 50333, 0, 3,
                                                                       48128, 31088, 48338,
                                                                       19938, 20088, 33413,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 50648, 0, 3,
                                                                       48338, 31238, 48548,
                                                                       20088, 20238, 33638,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 50963, 0, 3,
                                                                       48548, 31388, 48758,
                                                                       20238, 20388, 33863,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 51278, 0, 3,
                                                                       48968, 31988, 49178,
                                                                       20688, 20838, 34538,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 51593, 0, 3,
                                                                       49178, 32138, 49388,
                                                                       20838, 20988, 34763,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 51908, 0, 3,
                                                                       49388, 32288, 49598,
                                                                       20988, 21138, 34988,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 52223, 0, 3,
                                                                       49598, 32438, 49808,
                                                                       21138, 21288, 35213,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 52538, 0, 3,
                                                                       50018, 33188, 50333,
                                                                       21588, 21798, 36068,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 52979, 0, 3,
                                                                       50333, 33413, 50648,
                                                                       21798, 22008, 36383,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 53420, 0, 3,
                                                                       50648, 33638, 50963,
                                                                       22008, 22218, 36698,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 53861, 0, 3,
                                                                       51278, 34538, 51593,
                                                                       22638, 22848, 37643,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 54302, 0, 3,
                                                                       51593, 34763, 51908,
                                                                       22848, 23058, 37958,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsh_three_center_electron_repulsion_0(buffer, 54743, 0, 3,
                                                                       51908, 34988, 52223,
                                                                       23058, 23268, 38273,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 55184, 0, 3,
                                                                       52538, 36068, 52979,
                                                                       23688, 23968, 39428,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 55772, 0, 3,
                                                                       52979, 36383, 53420,
                                                                       23968, 24248, 39848,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 56360, 0, 3,
                                                                       53861, 37643, 54302,
                                                                       24808, 25088, 41108,
                                                                       ncols, gamma, p, q);

                    compute_prim_ish_three_center_electron_repulsion_0(buffer, 56948, 0, 3,
                                                                       54302, 37958, 54743,
                                                                       25088, 25368, 41528,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 57536, 0, 3,
                                                                       55184, 39428, 55772,
                                                                       25928, 26288, 43028,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksh_three_center_electron_repulsion_0(buffer, 58292, 0, 3,
                                                                       56360, 41108, 56948,
                                                                       27008, 27368, 44648,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59048, 3, 28088,
                                                                       28103, 45188, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59076, 3, 28103,
                                                                       28118, 45209, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59104, 3, 28118,
                                                                       28133, 45230, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59132, 3, 28133,
                                                                       28148, 45251, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59160, 3, 28148,
                                                                       28163, 45272, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59188, 3, 28163,
                                                                       28178, 45293, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59216, 3, 28178,
                                                                       28193, 45314, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59244, 3, 28193,
                                                                       28208, 45335, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59272, 3, 28238,
                                                                       28253, 45356, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59300, 3, 28253,
                                                                       28268, 45377, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59328, 3, 28268,
                                                                       28283, 45398, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59356, 3, 28283,
                                                                       28298, 45419, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59384, 3, 28298,
                                                                       28313, 45440, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59412, 3, 28313,
                                                                       28328, 45461, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59440, 3, 28328,
                                                                       28343, 45482, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 59468, 3, 28343,
                                                                       28358, 45503, ncols,
                                                                       gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 59496, 0, 3,
                                                                       59048, 45188, 59076,
                                                                       28388, 28433, 45524,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 59580, 0, 3,
                                                                       59076, 45209, 59104,
                                                                       28433, 28478, 45587,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 59664, 0, 3,
                                                                       59104, 45230, 59132,
                                                                       28478, 28523, 45650,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 59748, 0, 3,
                                                                       59132, 45251, 59160,
                                                                       28523, 28568, 45713,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 59832, 0, 3,
                                                                       59160, 45272, 59188,
                                                                       28568, 28613, 45776,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 59916, 0, 3,
                                                                       59188, 45293, 59216,
                                                                       28613, 28658, 45839,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 60000, 0, 3,
                                                                       59216, 45314, 59244,
                                                                       28658, 28703, 45902,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 60084, 0, 3,
                                                                       59272, 45356, 59300,
                                                                       28793, 28838, 45965,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 60168, 0, 3,
                                                                       59300, 45377, 59328,
                                                                       28838, 28883, 46028,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 60252, 0, 3,
                                                                       59328, 45398, 59356,
                                                                       28883, 28928, 46091,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 60336, 0, 3,
                                                                       59356, 45419, 59384,
                                                                       28928, 28973, 46154,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 60420, 0, 3,
                                                                       59384, 45440, 59412,
                                                                       28973, 29018, 46217,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 60504, 0, 3,
                                                                       59412, 45461, 59440,
                                                                       29018, 29063, 46280,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 60588, 0, 3,
                                                                       59440, 45482, 59468,
                                                                       29063, 29108, 46343,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 60672, 0, 3,
                                                                       59496, 45524, 59580,
                                                                       29198, 29288, 46406,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 60840, 0, 3,
                                                                       59580, 45587, 59664,
                                                                       29288, 29378, 46532,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 61008, 0, 3,
                                                                       59664, 45650, 59748,
                                                                       29378, 29468, 46658,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 61176, 0, 3,
                                                                       59748, 45713, 59832,
                                                                       29468, 29558, 46784,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 61344, 0, 3,
                                                                       59832, 45776, 59916,
                                                                       29558, 29648, 46910,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 61512, 0, 3,
                                                                       59916, 45839, 60000,
                                                                       29648, 29738, 47036,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 61680, 0, 3,
                                                                       60084, 45965, 60168,
                                                                       29918, 30008, 47162,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 61848, 0, 3,
                                                                       60168, 46028, 60252,
                                                                       30008, 30098, 47288,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 62016, 0, 3,
                                                                       60252, 46091, 60336,
                                                                       30098, 30188, 47414,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 62184, 0, 3,
                                                                       60336, 46154, 60420,
                                                                       30188, 30278, 47540,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 62352, 0, 3,
                                                                       60420, 46217, 60504,
                                                                       30278, 30368, 47666,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 62520, 0, 3,
                                                                       60504, 46280, 60588,
                                                                       30368, 30458, 47792,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 62688, 0, 3,
                                                                       60672, 46406, 60840,
                                                                       30638, 30788, 47918,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 62968, 0, 3,
                                                                       60840, 46532, 61008,
                                                                       30788, 30938, 48128,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 63248, 0, 3,
                                                                       61008, 46658, 61176,
                                                                       30938, 31088, 48338,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 63528, 0, 3,
                                                                       61176, 46784, 61344,
                                                                       31088, 31238, 48548,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 63808, 0, 3,
                                                                       61344, 46910, 61512,
                                                                       31238, 31388, 48758,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 64088, 0, 3,
                                                                       61680, 47162, 61848,
                                                                       31688, 31838, 48968,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 64368, 0, 3,
                                                                       61848, 47288, 62016,
                                                                       31838, 31988, 49178,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 64648, 0, 3,
                                                                       62016, 47414, 62184,
                                                                       31988, 32138, 49388,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 64928, 0, 3,
                                                                       62184, 47540, 62352,
                                                                       32138, 32288, 49598,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 65208, 0, 3,
                                                                       62352, 47666, 62520,
                                                                       32288, 32438, 49808,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 65488, 0, 3,
                                                                       62688, 47918, 62968,
                                                                       32738, 32963, 50018,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 65908, 0, 3,
                                                                       62968, 48128, 63248,
                                                                       32963, 33188, 50333,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 66328, 0, 3,
                                                                       63248, 48338, 63528,
                                                                       33188, 33413, 50648,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 66748, 0, 3,
                                                                       63528, 48548, 63808,
                                                                       33413, 33638, 50963,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 67168, 0, 3,
                                                                       64088, 48968, 64368,
                                                                       34088, 34313, 51278,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 67588, 0, 3,
                                                                       64368, 49178, 64648,
                                                                       34313, 34538, 51593,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 68008, 0, 3,
                                                                       64648, 49388, 64928,
                                                                       34538, 34763, 51908,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 68428, 0, 3,
                                                                       64928, 49598, 65208,
                                                                       34763, 34988, 52223,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 68848, 0, 3,
                                                                       65488, 50018, 65908,
                                                                       35438, 35753, 52538,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 69436, 0, 3,
                                                                       65908, 50333, 66328,
                                                                       35753, 36068, 52979,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 70024, 0, 3,
                                                                       66328, 50648, 66748,
                                                                       36068, 36383, 53420,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 70612, 0, 3,
                                                                       67168, 51278, 67588,
                                                                       37013, 37328, 53861,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 71200, 0, 3,
                                                                       67588, 51593, 68008,
                                                                       37328, 37643, 54302,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsi_three_center_electron_repulsion_0(buffer, 71788, 0, 3,
                                                                       68008, 51908, 68428,
                                                                       37643, 37958, 54743,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 72376, 0, 3,
                                                                       68848, 52538, 69436,
                                                                       38588, 39008, 55184,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 73160, 0, 3,
                                                                       69436, 52979, 70024,
                                                                       39008, 39428, 55772,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 73944, 0, 3,
                                                                       70612, 53861, 71200,
                                                                       40268, 40688, 56360,
                                                                       ncols, gamma, p, q);

                    compute_prim_isi_three_center_electron_repulsion_0(buffer, 74728, 0, 3,
                                                                       71200, 54302, 71788,
                                                                       40688, 41108, 56948,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 75512, 0, 3,
                                                                       72376, 55184, 73160,
                                                                       41948, 42488, 57536,
                                                                       ncols, gamma, p, q);

                    compute_prim_ksi_three_center_electron_repulsion_0(buffer, 76520, 0, 3,
                                                                       73944, 56360, 74728,
                                                                       43568, 44108, 58292,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_f_x(buffer, 77528, 60672, 65488, 1, 28, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 77808, 60672, 65488, 1, 28, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 78088, 60672, 65488, 1, 28, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 78368, 61680, 67168, 1, 28, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 78648, 61680, 67168, 1, 28, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 78928, 61680, 67168, 1, 28, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 79208, 62688, 68848, 1, 28, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 79628, 62688, 68848, 1, 28, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 80048, 62688, 68848, 1, 28, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 80468, 64088, 70612, 1, 28, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 80888, 64088, 70612, 1, 28, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 81308, 64088, 70612, 1, 28, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 81728, 65488, 72376, 1, 28, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 82316, 65488, 72376, 1, 28, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 82904, 65488, 72376, 1, 28, ncols, alpha);

                    simdgeo::geom_h_x(buffer, 83492, 67168, 73944, 1, 28, ncols, alpha);

                    simdgeo::geom_h_y(buffer, 84080, 67168, 73944, 1, 28, ncols, alpha);

                    simdgeo::geom_h_z(buffer, 84668, 67168, 73944, 1, 28, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 85256, 68848, 75512, 1, 28, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 86040, 68848, 75512, 1, 28, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 86824, 68848, 75512, 1, 28, ncols, alpha);

                    simdgeo::geom_i_x(buffer, 87608, 70612, 76520, 1, 28, ncols, alpha);

                    simdgeo::geom_i_y(buffer, 88392, 70612, 76520, 1, 28, ncols, alpha);

                    simdgeo::geom_i_z(buffer, 89176, 70612, 76520, 1, 28, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 89960, 77528, 280, ncols);

                    simdfunc::contract_primitives(buffer, 90370, 77808, 280, ncols);

                    simdfunc::contract_primitives(buffer, 90780, 78088, 280, ncols);

                    simdfunc::contract_primitives(buffer, 91190, 62688, 280, ncols);

                    simdfunc::contract_primitives(buffer, 91600, 78368, 280, ncols);

                    simdfunc::contract_primitives(buffer, 92010, 78648, 280, ncols);

                    simdfunc::contract_primitives(buffer, 92420, 78928, 280, ncols);

                    simdfunc::contract_primitives(buffer, 92830, 64088, 280, ncols);

                    simdfunc::contract_primitives(buffer, 93240, 79208, 420, ncols);

                    simdfunc::contract_primitives(buffer, 93855, 79628, 420, ncols);

                    simdfunc::contract_primitives(buffer, 94470, 80048, 420, ncols);

                    simdfunc::contract_primitives(buffer, 95085, 65488, 420, ncols);

                    simdfunc::contract_primitives(buffer, 95700, 80468, 420, ncols);

                    simdfunc::contract_primitives(buffer, 96315, 80888, 420, ncols);

                    simdfunc::contract_primitives(buffer, 96930, 81308, 420, ncols);

                    simdfunc::contract_primitives(buffer, 97545, 67168, 420, ncols);

                    simdfunc::contract_primitives(buffer, 98160, 81728, 588, ncols);

                    simdfunc::contract_primitives(buffer, 99021, 82316, 588, ncols);

                    simdfunc::contract_primitives(buffer, 99882, 82904, 588, ncols);

                    simdfunc::contract_primitives(buffer, 100743, 68848, 588, ncols);

                    simdfunc::contract_primitives(buffer, 101604, 83492, 588, ncols);

                    simdfunc::contract_primitives(buffer, 102465, 84080, 588, ncols);

                    simdfunc::contract_primitives(buffer, 103326, 84668, 588, ncols);

                    simdfunc::contract_primitives(buffer, 104187, 70612, 588, ncols);

                    simdfunc::contract_primitives(buffer, 105048, 85256, 784, ncols);

                    simdfunc::contract_primitives(buffer, 106196, 86040, 784, ncols);

                    simdfunc::contract_primitives(buffer, 107344, 86824, 784, ncols);

                    simdfunc::contract_primitives(buffer, 108492, 87608, 784, ncols);

                    simdfunc::contract_primitives(buffer, 109640, 88392, 784, ncols);

                    simdfunc::contract_primitives(buffer, 110788, 89176, 784, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 90240, 89960, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 90650, 90370, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 91060, 90780, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 91470, 91190, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 91880, 91600, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 92290, 92010, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 92700, 92420, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 93110, 92830, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 93660, 93240, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 94275, 93855, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 94890, 94470, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 95505, 95085, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 96120, 95700, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 96735, 96315, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 97350, 96930, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 97965, 97545, 15, 1, nmax);

        simdtrf::transform_i_inner(buffer, 98748, 98160, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 99609, 99021, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 100470, 99882, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 101331, 100743, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 102192, 101604, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 103053, 102465, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 103914, 103326, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 104775, 104187, 21, 1, nmax);

        simdtrf::transform_i_inner(buffer, 105832, 105048, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 106980, 106196, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 108128, 107344, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 109276, 108492, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 110424, 109640, 28, 1, nmax);

        simdtrf::transform_i_inner(buffer, 111572, 110788, 28, 1, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 111936, 90240, 91470,
                                                       93660, 13, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 112326, 90650, 91470,
                                                       94275, 13, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 112716, 91060, 91470,
                                                       94890, 13, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 113106, 91470, 95505, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 113496, 91880, 93110,
                                                       96120, 13, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 113886, 92290, 93110,
                                                       96735, 13, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 114276, 92700, 93110,
                                                       97350, 13, nmax);

        simdtrf::compute_hrr_fp_out_of_first(buffer, coordinates, 114666, 93110, 97965, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 115056, 93660, 95505,
                                                       98748, 13, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 115641, 94275, 95505,
                                                       99609, 13, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 116226, 94890, 95505,
                                                       100470, 13, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 116811, 95505, 101331, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gp_out_of_first(buffer, coordinates, 117396, 96120, 97965,
                                                       102192, 13, nmax);

        simdtrf::compute_hrr_geom_100y_gp_out_of_first(buffer, coordinates, 117981, 96735, 97965,
                                                       103053, 13, nmax);

        simdtrf::compute_hrr_geom_100z_gp_out_of_first(buffer, coordinates, 118566, 97350, 97965,
                                                       103914, 13, nmax);

        simdtrf::compute_hrr_gp_out_of_first(buffer, coordinates, 119151, 97965, 104775, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 119736, 98748,
                                                       101331, 105832, 13, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 120555, 99609,
                                                       101331, 106980, 13, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 121374, 100470,
                                                       101331, 108128, 13, nmax);

        simdtrf::compute_hrr_geom_100x_hp_out_of_first(buffer, coordinates, 122193, 102192,
                                                       104775, 109276, 13, nmax);

        simdtrf::compute_hrr_geom_100y_hp_out_of_first(buffer, coordinates, 123012, 103053,
                                                       104775, 110424, 13, nmax);

        simdtrf::compute_hrr_geom_100z_hp_out_of_first(buffer, coordinates, 123831, 103914,
                                                       104775, 111572, 13, nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 124650, 111936,
                                                       113106, 115056, 13, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 125430, 112326,
                                                       113106, 115641, 13, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 126210, 112716,
                                                       113106, 116226, 13, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 126990, 113106, 116811, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_fd_out_of_first(buffer, coordinates, 127770, 113496,
                                                       114666, 117396, 13, nmax);

        simdtrf::compute_hrr_geom_100y_fd_out_of_first(buffer, coordinates, 128550, 113886,
                                                       114666, 117981, 13, nmax);

        simdtrf::compute_hrr_geom_100z_fd_out_of_first(buffer, coordinates, 129330, 114276,
                                                       114666, 118566, 13, nmax);

        simdtrf::compute_hrr_fd_out_of_first(buffer, coordinates, 130110, 114666, 119151, 13,
                                             nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 130890, 115056,
                                                       116811, 119736, 13, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 132060, 115641,
                                                       116811, 120555, 13, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 133230, 116226,
                                                       116811, 121374, 13, nmax);

        simdtrf::compute_hrr_geom_100x_gd_out_of_first(buffer, coordinates, 134400, 117396,
                                                       119151, 122193, 13, nmax);

        simdtrf::compute_hrr_geom_100y_gd_out_of_first(buffer, coordinates, 135570, 117981,
                                                       119151, 123012, 13, nmax);

        simdtrf::compute_hrr_geom_100z_gd_out_of_first(buffer, coordinates, 136740, 118566,
                                                       119151, 123831, 13, nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 137910, 124650,
                                                       126990, 130890, 13, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 139210, 125430,
                                                       126990, 132060, 13, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 140510, 126210,
                                                       126990, 133230, 13, nmax);

        simdtrf::compute_hrr_geom_100x_ff_out_of_first(buffer, coordinates, 141810, 127770,
                                                       130110, 134400, 13, nmax);

        simdtrf::compute_hrr_geom_100y_ff_out_of_first(buffer, coordinates, 143110, 128550,
                                                       130110, 135570, 13, nmax);

        simdtrf::compute_hrr_geom_100z_ff_out_of_first(buffer, coordinates, 144410, 129330,
                                                       130110, 136740, 13, nmax);

        simdtrf::transform_f_inner(buffer, 145710, 141810, 10, 13, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 145710, 91, nmax);

        simdtrf::transform_f_inner(buffer, 145710, 143110, 10, 13, nmax);

        simdtrf::transform_f_outer(values + 637 * nvalues + n * npairs, nvalues, buffer, 145710,
                                   91, nmax);

        simdtrf::transform_f_inner(buffer, 145710, 144410, 10, 13, nmax);

        simdtrf::transform_f_outer(values + 1274 * nvalues + n * npairs, nvalues, buffer, 145710,
                                   91, nmax);

        simdtrf::transform_f_inner(buffer, 145710, 137910, 10, 13, nmax);

        simdtrf::transform_f_outer(values + 1911 * nvalues + n * npairs, nvalues, buffer, 145710,
                                   91, nmax);

        simdtrf::transform_f_inner(buffer, 145710, 139210, 10, 13, nmax);

        simdtrf::transform_f_outer(values + 2548 * nvalues + n * npairs, nvalues, buffer, 145710,
                                   91, nmax);

        simdtrf::transform_f_inner(buffer, 145710, 140510, 10, 13, nmax);

        simdtrf::transform_f_outer(values + 3185 * nvalues + n * npairs, nvalues, buffer, 145710,
                                   91, nmax);
    }

    for (size_t m = 0; m < 3822; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
