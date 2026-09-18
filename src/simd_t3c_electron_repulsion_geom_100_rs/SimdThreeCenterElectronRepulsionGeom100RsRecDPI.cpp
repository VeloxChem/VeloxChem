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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecDPI.hpp"

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
#include "SimdTransferGeom100XDP.hpp"
#include "SimdTransferGeom100YDP.hpp"
#include "SimdTransferGeom100ZDP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_dpi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_dpi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 26962, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 1170 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 26962, 20896, 4298, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 568, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 571, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 574, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 577, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 580, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 583, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 586, 3, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 589, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 592, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 595, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 598, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 601, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 604, 3, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 607, 3, 25, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 610, 3, 26, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 613, 3, 27, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 616, 3, 28, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 619, 3, 29, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 622, 3, 9, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 631, 3, 10, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 640, 3, 11, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 649, 3, 12, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 658, 3, 13, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 667, 3, 14, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 676, 3, 15, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 685, 3, 16, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 694, 3, 21, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 703, 3, 22, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 712, 3, 23, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 721, 3, 24, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 730, 3, 25, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 739, 3, 26, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 748, 3, 27, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 757, 3, 28, 87,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 766, 3, 36, 102,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 784, 3, 39, 108,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 802, 3, 42, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 820, 3, 45, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 838, 3, 48, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 856, 3, 51, 132,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 874, 3, 54, 138,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 892, 3, 66, 156,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 910, 3, 69, 162,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 928, 3, 72, 168,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 946, 3, 75, 174,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 964, 3, 78, 180,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 982, 3, 81, 186,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1000, 3, 84, 192,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1018, 3, 102, 218,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1048, 3, 108, 228,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1078, 3, 114, 238,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1108, 3, 120, 248,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1138, 3, 126, 258,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1168, 3, 132, 268,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1198, 3, 156, 298,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1228, 3, 162, 308,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1258, 3, 168, 318,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1288, 3, 174, 328,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1318, 3, 180, 338,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1348, 3, 186, 348,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1378, 3, 218, 388,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1423, 3, 228, 403,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1468, 3, 238, 418,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1513, 3, 248, 433,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1558, 3, 258, 448,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1603, 3, 298, 493,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1648, 3, 308, 508,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1693, 3, 318, 523,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1738, 3, 328, 538,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1783, 3, 338, 553,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1828, 3, 7, 8,
                                                                       568, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1834, 3, 8, 9,
                                                                       571, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1840, 3, 9, 10,
                                                                       574, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1846, 3, 10, 11,
                                                                       577, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1852, 3, 11, 12,
                                                                       580, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1858, 3, 12, 13,
                                                                       583, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1864, 3, 13, 14,
                                                                       586, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1870, 3, 14, 15,
                                                                       589, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1876, 3, 15, 16,
                                                                       592, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1882, 3, 19, 20,
                                                                       595, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1888, 3, 20, 21,
                                                                       598, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1894, 3, 21, 22,
                                                                       601, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1900, 3, 22, 23,
                                                                       604, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1906, 3, 23, 24,
                                                                       607, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1912, 3, 24, 25,
                                                                       610, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1918, 3, 25, 26,
                                                                       613, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1924, 3, 26, 27,
                                                                       616, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1930, 3, 27, 28,
                                                                       619, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1936, 0, 3, 1828,
                                                                       568, 1834, 622, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1954, 0, 3, 1834,
                                                                       571, 1840, 631, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1972, 0, 3, 1840,
                                                                       574, 1846, 640, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1990, 0, 3, 1846,
                                                                       577, 1852, 649, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2008, 0, 3, 1852,
                                                                       580, 1858, 658, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2026, 0, 3, 1858,
                                                                       583, 1864, 667, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2044, 0, 3, 1864,
                                                                       586, 1870, 676, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2062, 0, 3, 1870,
                                                                       589, 1876, 685, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2080, 0, 3, 1882,
                                                                       595, 1888, 694, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2098, 0, 3, 1888,
                                                                       598, 1894, 703, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2116, 0, 3, 1894,
                                                                       601, 1900, 712, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2134, 0, 3, 1900,
                                                                       604, 1906, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2152, 0, 3, 1906,
                                                                       607, 1912, 730, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2170, 0, 3, 1912,
                                                                       610, 1918, 739, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2188, 0, 3, 1918,
                                                                       613, 1924, 748, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2206, 0, 3, 1924,
                                                                       616, 1930, 757, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2224, 0, 3, 1936,
                                                                       622, 1954, 90, 96, 766,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2260, 0, 3, 1954,
                                                                       631, 1972, 96, 102, 784,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2296, 0, 3, 1972,
                                                                       640, 1990, 102, 108, 802,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2332, 0, 3, 1990,
                                                                       649, 2008, 108, 114, 820,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2368, 0, 3, 2008,
                                                                       658, 2026, 114, 120, 838,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2404, 0, 3, 2026,
                                                                       667, 2044, 120, 126, 856,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2440, 0, 3, 2044,
                                                                       676, 2062, 126, 132, 874,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2476, 0, 3, 2080,
                                                                       694, 2098, 144, 150, 892,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2512, 0, 3, 2098,
                                                                       703, 2116, 150, 156, 910,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2548, 0, 3, 2116,
                                                                       712, 2134, 156, 162, 928,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2584, 0, 3, 2134,
                                                                       721, 2152, 162, 168, 946,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2620, 0, 3, 2152,
                                                                       730, 2170, 168, 174, 964,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2656, 0, 3, 2170,
                                                                       739, 2188, 174, 180, 982,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2692, 0, 3, 2188,
                                                                       748, 2206, 180, 186, 1000,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2728, 0, 3, 2224,
                                                                       766, 2260, 198, 208, 1018,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2788, 0, 3, 2260,
                                                                       784, 2296, 208, 218, 1048,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2848, 0, 3, 2296,
                                                                       802, 2332, 218, 228, 1078,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2908, 0, 3, 2332,
                                                                       820, 2368, 228, 238, 1108,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2968, 0, 3, 2368,
                                                                       838, 2404, 238, 248, 1138,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3028, 0, 3, 2404,
                                                                       856, 2440, 248, 258, 1168,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3088, 0, 3, 2476,
                                                                       892, 2512, 278, 288, 1198,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3148, 0, 3, 2512,
                                                                       910, 2548, 288, 298, 1228,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3208, 0, 3, 2548,
                                                                       928, 2584, 298, 308, 1258,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3268, 0, 3, 2584,
                                                                       946, 2620, 308, 318, 1288,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3328, 0, 3, 2620,
                                                                       964, 2656, 318, 328, 1318,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3388, 0, 3, 2656,
                                                                       982, 2692, 328, 338, 1348,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3448, 0, 3, 2728,
                                                                       1018, 2788, 358, 373,
                                                                       1378, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3538, 0, 3, 2788,
                                                                       1048, 2848, 373, 388,
                                                                       1423, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3628, 0, 3, 2848,
                                                                       1078, 2908, 388, 403,
                                                                       1468, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3718, 0, 3, 2908,
                                                                       1108, 2968, 403, 418,
                                                                       1513, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3808, 0, 3, 2968,
                                                                       1138, 3028, 418, 433,
                                                                       1558, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3898, 0, 3, 3088,
                                                                       1198, 3148, 463, 478,
                                                                       1603, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3988, 0, 3, 3148,
                                                                       1228, 3208, 478, 493,
                                                                       1648, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4078, 0, 3, 3208,
                                                                       1258, 3268, 493, 508,
                                                                       1693, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4168, 0, 3, 3268,
                                                                       1288, 3328, 508, 523,
                                                                       1738, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 4258, 0, 3, 3328,
                                                                       1318, 3388, 523, 538,
                                                                       1783, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4348, 3, 568, 571,
                                                                       1840, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4358, 3, 571, 574,
                                                                       1846, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4368, 3, 574, 577,
                                                                       1852, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4378, 3, 577, 580,
                                                                       1858, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4388, 3, 580, 583,
                                                                       1864, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4398, 3, 583, 586,
                                                                       1870, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4408, 3, 586, 589,
                                                                       1876, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4418, 3, 595, 598,
                                                                       1894, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4428, 3, 598, 601,
                                                                       1900, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4438, 3, 601, 604,
                                                                       1906, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4448, 3, 604, 607,
                                                                       1912, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4458, 3, 607, 610,
                                                                       1918, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4468, 3, 610, 613,
                                                                       1924, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4478, 3, 613, 616,
                                                                       1930, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4488, 0, 3, 4348,
                                                                       1840, 4358, 622, 631,
                                                                       1972, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4518, 0, 3, 4358,
                                                                       1846, 4368, 631, 640,
                                                                       1990, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4548, 0, 3, 4368,
                                                                       1852, 4378, 640, 649,
                                                                       2008, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4578, 0, 3, 4378,
                                                                       1858, 4388, 649, 658,
                                                                       2026, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4608, 0, 3, 4388,
                                                                       1864, 4398, 658, 667,
                                                                       2044, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4638, 0, 3, 4398,
                                                                       1870, 4408, 667, 676,
                                                                       2062, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4668, 0, 3, 4418,
                                                                       1894, 4428, 694, 703,
                                                                       2116, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4698, 0, 3, 4428,
                                                                       1900, 4438, 703, 712,
                                                                       2134, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4728, 0, 3, 4438,
                                                                       1906, 4448, 712, 721,
                                                                       2152, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4758, 0, 3, 4448,
                                                                       1912, 4458, 721, 730,
                                                                       2170, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4788, 0, 3, 4458,
                                                                       1918, 4468, 730, 739,
                                                                       2188, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4818, 0, 3, 4468,
                                                                       1924, 4478, 739, 748,
                                                                       2206, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4848, 0, 3, 4488,
                                                                       1972, 4518, 766, 784,
                                                                       2296, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4908, 0, 3, 4518,
                                                                       1990, 4548, 784, 802,
                                                                       2332, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4968, 0, 3, 4548,
                                                                       2008, 4578, 802, 820,
                                                                       2368, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5028, 0, 3, 4578,
                                                                       2026, 4608, 820, 838,
                                                                       2404, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5088, 0, 3, 4608,
                                                                       2044, 4638, 838, 856,
                                                                       2440, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5148, 0, 3, 4668,
                                                                       2116, 4698, 892, 910,
                                                                       2548, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5208, 0, 3, 4698,
                                                                       2134, 4728, 910, 928,
                                                                       2584, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5268, 0, 3, 4728,
                                                                       2152, 4758, 928, 946,
                                                                       2620, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5328, 0, 3, 4758,
                                                                       2170, 4788, 946, 964,
                                                                       2656, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5388, 0, 3, 4788,
                                                                       2188, 4818, 964, 982,
                                                                       2692, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5448, 0, 3, 4848,
                                                                       2296, 4908, 1018, 1048,
                                                                       2848, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5548, 0, 3, 4908,
                                                                       2332, 4968, 1048, 1078,
                                                                       2908, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5648, 0, 3, 4968,
                                                                       2368, 5028, 1078, 1108,
                                                                       2968, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5748, 0, 3, 5028,
                                                                       2404, 5088, 1108, 1138,
                                                                       3028, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5848, 0, 3, 5148,
                                                                       2548, 5208, 1198, 1228,
                                                                       3208, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 5948, 0, 3, 5208,
                                                                       2584, 5268, 1228, 1258,
                                                                       3268, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6048, 0, 3, 5268,
                                                                       2620, 5328, 1258, 1288,
                                                                       3328, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6148, 0, 3, 5328,
                                                                       2656, 5388, 1288, 1318,
                                                                       3388, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 6248, 0, 3, 5448,
                                                                       2848, 5548, 1378, 1423,
                                                                       3628, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 6398, 0, 3, 5548,
                                                                       2908, 5648, 1423, 1468,
                                                                       3718, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 6548, 0, 3, 5648,
                                                                       2968, 5748, 1468, 1513,
                                                                       3808, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 6698, 0, 3, 5848,
                                                                       3208, 5948, 1603, 1648,
                                                                       4078, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 6848, 0, 3, 5948,
                                                                       3268, 6048, 1648, 1693,
                                                                       4168, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 6998, 0, 3, 6048,
                                                                       3328, 6148, 1693, 1738,
                                                                       4258, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7148, 3, 1828,
                                                                       1834, 4348, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7163, 3, 1834,
                                                                       1840, 4358, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7178, 3, 1840,
                                                                       1846, 4368, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7193, 3, 1846,
                                                                       1852, 4378, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7208, 3, 1852,
                                                                       1858, 4388, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7223, 3, 1858,
                                                                       1864, 4398, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7238, 3, 1864,
                                                                       1870, 4408, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7253, 3, 1882,
                                                                       1888, 4418, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7268, 3, 1888,
                                                                       1894, 4428, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7283, 3, 1894,
                                                                       1900, 4438, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7298, 3, 1900,
                                                                       1906, 4448, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7313, 3, 1906,
                                                                       1912, 4458, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7328, 3, 1912,
                                                                       1918, 4468, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7343, 3, 1918,
                                                                       1924, 4478, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7358, 0, 3, 7148,
                                                                       4348, 7163, 1936, 1954,
                                                                       4488, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7403, 0, 3, 7163,
                                                                       4358, 7178, 1954, 1972,
                                                                       4518, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7448, 0, 3, 7178,
                                                                       4368, 7193, 1972, 1990,
                                                                       4548, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7493, 0, 3, 7193,
                                                                       4378, 7208, 1990, 2008,
                                                                       4578, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7538, 0, 3, 7208,
                                                                       4388, 7223, 2008, 2026,
                                                                       4608, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7583, 0, 3, 7223,
                                                                       4398, 7238, 2026, 2044,
                                                                       4638, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7628, 0, 3, 7253,
                                                                       4418, 7268, 2080, 2098,
                                                                       4668, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7673, 0, 3, 7268,
                                                                       4428, 7283, 2098, 2116,
                                                                       4698, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7718, 0, 3, 7283,
                                                                       4438, 7298, 2116, 2134,
                                                                       4728, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7763, 0, 3, 7298,
                                                                       4448, 7313, 2134, 2152,
                                                                       4758, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7808, 0, 3, 7313,
                                                                       4458, 7328, 2152, 2170,
                                                                       4788, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7853, 0, 3, 7328,
                                                                       4468, 7343, 2170, 2188,
                                                                       4818, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 7898, 0, 3, 7358,
                                                                       4488, 7403, 2224, 2260,
                                                                       4848, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 7988, 0, 3, 7403,
                                                                       4518, 7448, 2260, 2296,
                                                                       4908, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8078, 0, 3, 7448,
                                                                       4548, 7493, 2296, 2332,
                                                                       4968, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8168, 0, 3, 7493,
                                                                       4578, 7538, 2332, 2368,
                                                                       5028, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8258, 0, 3, 7538,
                                                                       4608, 7583, 2368, 2404,
                                                                       5088, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8348, 0, 3, 7628,
                                                                       4668, 7673, 2476, 2512,
                                                                       5148, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8438, 0, 3, 7673,
                                                                       4698, 7718, 2512, 2548,
                                                                       5208, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8528, 0, 3, 7718,
                                                                       4728, 7763, 2548, 2584,
                                                                       5268, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8618, 0, 3, 7763,
                                                                       4758, 7808, 2584, 2620,
                                                                       5328, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 8708, 0, 3, 7808,
                                                                       4788, 7853, 2620, 2656,
                                                                       5388, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 8798, 0, 3, 7898,
                                                                       4848, 7988, 2728, 2788,
                                                                       5448, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 8948, 0, 3, 7988,
                                                                       4908, 8078, 2788, 2848,
                                                                       5548, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 9098, 0, 3, 8078,
                                                                       4968, 8168, 2848, 2908,
                                                                       5648, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 9248, 0, 3, 8168,
                                                                       5028, 8258, 2908, 2968,
                                                                       5748, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 9398, 0, 3, 8348,
                                                                       5148, 8438, 3088, 3148,
                                                                       5848, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 9548, 0, 3, 8438,
                                                                       5208, 8528, 3148, 3208,
                                                                       5948, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 9698, 0, 3, 8528,
                                                                       5268, 8618, 3208, 3268,
                                                                       6048, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 9848, 0, 3, 8618,
                                                                       5328, 8708, 3268, 3328,
                                                                       6148, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 9998, 0, 3, 8798,
                                                                       5448, 8948, 3448, 3538,
                                                                       6248, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 10223, 0, 3, 8948,
                                                                       5548, 9098, 3538, 3628,
                                                                       6398, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 10448, 0, 3, 9098,
                                                                       5648, 9248, 3628, 3718,
                                                                       6548, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 10673, 0, 3, 9398,
                                                                       5848, 9548, 3898, 3988,
                                                                       6698, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 10898, 0, 3, 9548,
                                                                       5948, 9698, 3988, 4078,
                                                                       6848, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 11123, 0, 3, 9698,
                                                                       6048, 9848, 4078, 4168,
                                                                       6998, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11348, 3, 4348,
                                                                       4358, 7178, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11369, 3, 4358,
                                                                       4368, 7193, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11390, 3, 4368,
                                                                       4378, 7208, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11411, 3, 4378,
                                                                       4388, 7223, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11432, 3, 4388,
                                                                       4398, 7238, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11453, 3, 4418,
                                                                       4428, 7283, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11474, 3, 4428,
                                                                       4438, 7298, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11495, 3, 4438,
                                                                       4448, 7313, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11516, 3, 4448,
                                                                       4458, 7328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 11537, 3, 4458,
                                                                       4468, 7343, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11558, 0, 3,
                                                                       11348, 7178, 11369, 4488,
                                                                       4518, 7448, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11621, 0, 3,
                                                                       11369, 7193, 11390, 4518,
                                                                       4548, 7493, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11684, 0, 3,
                                                                       11390, 7208, 11411, 4548,
                                                                       4578, 7538, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11747, 0, 3,
                                                                       11411, 7223, 11432, 4578,
                                                                       4608, 7583, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11810, 0, 3,
                                                                       11453, 7283, 11474, 4668,
                                                                       4698, 7718, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11873, 0, 3,
                                                                       11474, 7298, 11495, 4698,
                                                                       4728, 7763, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11936, 0, 3,
                                                                       11495, 7313, 11516, 4728,
                                                                       4758, 7808, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 11999, 0, 3,
                                                                       11516, 7328, 11537, 4758,
                                                                       4788, 7853, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 12062, 0, 3,
                                                                       11558, 7448, 11621, 4848,
                                                                       4908, 8078, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 12188, 0, 3,
                                                                       11621, 7493, 11684, 4908,
                                                                       4968, 8168, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 12314, 0, 3,
                                                                       11684, 7538, 11747, 4968,
                                                                       5028, 8258, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 12440, 0, 3,
                                                                       11810, 7718, 11873, 5148,
                                                                       5208, 8528, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 12566, 0, 3,
                                                                       11873, 7763, 11936, 5208,
                                                                       5268, 8618, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 12692, 0, 3,
                                                                       11936, 7808, 11999, 5268,
                                                                       5328, 8708, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 12818, 0, 3,
                                                                       12062, 8078, 12188, 5448,
                                                                       5548, 9098, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 13028, 0, 3,
                                                                       12188, 8168, 12314, 5548,
                                                                       5648, 9248, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 13238, 0, 3,
                                                                       12440, 8528, 12566, 5848,
                                                                       5948, 9698, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 13448, 0, 3,
                                                                       12566, 8618, 12692, 5948,
                                                                       6048, 9848, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 13658, 0, 3,
                                                                       12818, 9098, 13028, 6248,
                                                                       6398, 10448, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 13973, 0, 3,
                                                                       13238, 9698, 13448, 6698,
                                                                       6848, 11123, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 14288, 3, 7148,
                                                                       7163, 11348, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 14316, 3, 7163,
                                                                       7178, 11369, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 14344, 3, 7178,
                                                                       7193, 11390, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 14372, 3, 7193,
                                                                       7208, 11411, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 14400, 3, 7208,
                                                                       7223, 11432, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 14428, 3, 7253,
                                                                       7268, 11453, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 14456, 3, 7268,
                                                                       7283, 11474, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 14484, 3, 7283,
                                                                       7298, 11495, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 14512, 3, 7298,
                                                                       7313, 11516, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 14540, 3, 7313,
                                                                       7328, 11537, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 14568, 0, 3,
                                                                       14288, 11348, 14316, 7358,
                                                                       7403, 11558, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 14652, 0, 3,
                                                                       14316, 11369, 14344, 7403,
                                                                       7448, 11621, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 14736, 0, 3,
                                                                       14344, 11390, 14372, 7448,
                                                                       7493, 11684, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 14820, 0, 3,
                                                                       14372, 11411, 14400, 7493,
                                                                       7538, 11747, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 14904, 0, 3,
                                                                       14428, 11453, 14456, 7628,
                                                                       7673, 11810, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 14988, 0, 3,
                                                                       14456, 11474, 14484, 7673,
                                                                       7718, 11873, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 15072, 0, 3,
                                                                       14484, 11495, 14512, 7718,
                                                                       7763, 11936, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 15156, 0, 3,
                                                                       14512, 11516, 14540, 7763,
                                                                       7808, 11999, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 15240, 0, 3,
                                                                       14568, 11558, 14652, 7898,
                                                                       7988, 12062, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 15408, 0, 3,
                                                                       14652, 11621, 14736, 7988,
                                                                       8078, 12188, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 15576, 0, 3,
                                                                       14736, 11684, 14820, 8078,
                                                                       8168, 12314, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 15744, 0, 3,
                                                                       14904, 11810, 14988, 8348,
                                                                       8438, 12440, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 15912, 0, 3,
                                                                       14988, 11873, 15072, 8438,
                                                                       8528, 12566, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 16080, 0, 3,
                                                                       15072, 11936, 15156, 8528,
                                                                       8618, 12692, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 16248, 0, 3,
                                                                       15240, 12062, 15408, 8798,
                                                                       8948, 12818, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 16528, 0, 3,
                                                                       15408, 12188, 15576, 8948,
                                                                       9098, 13028, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 16808, 0, 3,
                                                                       15744, 12440, 15912, 9398,
                                                                       9548, 13238, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 17088, 0, 3,
                                                                       15912, 12566, 16080, 9548,
                                                                       9698, 13448, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 17368, 0, 3,
                                                                       16248, 12818, 16528, 9998,
                                                                       10223, 13658, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 17788, 0, 3,
                                                                       16808, 13238, 17088,
                                                                       10673, 10898, 13973,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_d_x(buffer, 18208, 14568, 16248, 1, 28, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 18376, 14568, 16248, 1, 28, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 18544, 14568, 16248, 1, 28, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 18712, 14904, 16808, 1, 28, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 18880, 14904, 16808, 1, 28, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 19048, 14904, 16808, 1, 28, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 19216, 15240, 17368, 1, 28, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 19496, 15240, 17368, 1, 28, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 19776, 15240, 17368, 1, 28, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 20056, 15744, 17788, 1, 28, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 20336, 15744, 17788, 1, 28, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 20616, 15744, 17788, 1, 28, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 20896, 18208, 168, ncols);

                    simdfunc::contract_primitives(buffer, 21142, 18376, 168, ncols);

                    simdfunc::contract_primitives(buffer, 21388, 18544, 168, ncols);

                    simdfunc::contract_primitives(buffer, 21634, 15240, 168, ncols);

                    simdfunc::contract_primitives(buffer, 21880, 18712, 168, ncols);

                    simdfunc::contract_primitives(buffer, 22126, 18880, 168, ncols);

                    simdfunc::contract_primitives(buffer, 22372, 19048, 168, ncols);

                    simdfunc::contract_primitives(buffer, 22618, 15744, 168, ncols);

                    simdfunc::contract_primitives(buffer, 22864, 19216, 280, ncols);

                    simdfunc::contract_primitives(buffer, 23274, 19496, 280, ncols);

                    simdfunc::contract_primitives(buffer, 23684, 19776, 280, ncols);

                    simdfunc::contract_primitives(buffer, 24094, 20056, 280, ncols);

                    simdfunc::contract_primitives(buffer, 24504, 20336, 280, ncols);

                    simdfunc::contract_primitives(buffer, 24914, 20616, 280, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 21064, 20896, 6, 1, nmax);

        simdtrf::transform_i_inner(buffer, 21310, 21142, 6, 1, nmax);

        simdtrf::transform_i_inner(buffer, 21556, 21388, 6, 1, nmax);

        simdtrf::transform_i_inner(buffer, 21802, 21634, 6, 1, nmax);

        simdtrf::transform_i_inner(buffer, 22048, 21880, 6, 1, nmax);

        simdtrf::transform_i_inner(buffer, 22294, 22126, 6, 1, nmax);

        simdtrf::transform_i_inner(buffer, 22540, 22372, 6, 1, nmax);

        simdtrf::transform_i_inner(buffer, 22786, 22618, 6, 1, nmax);

        simdtrf::transform_i_inner(buffer, 23144, 22864, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 23554, 23274, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 23964, 23684, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 24374, 24094, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 24784, 24504, 10, 1, nmax);

        simdtrf::transform_i_inner(buffer, 25194, 24914, 10, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 25324, 21064, 21802,
                                                       23144, 13, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 25558, 21310, 21802,
                                                       23554, 13, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 25792, 21556, 21802,
                                                       23964, 13, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 26026, 22048, 22786,
                                                       24374, 13, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 26260, 22294, 22786,
                                                       24784, 13, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 26494, 22540, 22786,
                                                       25194, 13, nmax);

        simdtrf::transform_p_inner(buffer, 26728, 26026, 6, 13, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 26728, 39, nmax);

        simdtrf::transform_p_inner(buffer, 26728, 26260, 6, 13, nmax);

        simdtrf::transform_d_outer(values + 195 * nvalues + n * npairs, nvalues, buffer, 26728,
                                   39, nmax);

        simdtrf::transform_p_inner(buffer, 26728, 26494, 6, 13, nmax);

        simdtrf::transform_d_outer(values + 390 * nvalues + n * npairs, nvalues, buffer, 26728,
                                   39, nmax);

        simdtrf::transform_p_inner(buffer, 26728, 25324, 6, 13, nmax);

        simdtrf::transform_d_outer(values + 585 * nvalues + n * npairs, nvalues, buffer, 26728,
                                   39, nmax);

        simdtrf::transform_p_inner(buffer, 26728, 25558, 6, 13, nmax);

        simdtrf::transform_d_outer(values + 780 * nvalues + n * npairs, nvalues, buffer, 26728,
                                   39, nmax);

        simdtrf::transform_p_inner(buffer, 26728, 25792, 6, 13, nmax);

        simdtrf::transform_d_outer(values + 975 * nvalues + n * npairs, nvalues, buffer, 26728,
                                   39, nmax);
    }

    for (size_t m = 0; m < 1170; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
