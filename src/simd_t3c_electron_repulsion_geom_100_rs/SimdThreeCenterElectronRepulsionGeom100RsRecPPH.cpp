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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecPPH.hpp"

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
#include "SimdGeometryP1.hpp"
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
#include "SimdTransferGeom100XPP.hpp"
#include "SimdTransferGeom100YPP.hpp"
#include "SimdTransferGeom100ZPP.hpp"
#include "SimdTransformH.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_pph_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_pph_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 9595, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 594 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 9595, 6982, 1854, dimensions);

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
                                                            4, 5, 6, 7, 8}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 15, 3, {1, 2, 3, 4,
                                                        5, 6, 7, 8}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 24, 0, 3, 7, 8,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 27, 0, 3, 8, 9,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 30, 0, 3, 9, 10,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 33, 0, 3, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 36, 0, 3, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 39, 0, 3, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 42, 0, 3, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 45, 0, 3, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 48, 0, 3, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 51, 0, 3, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 54, 0, 3, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 57, 0, 3, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 60, 0, 3, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 63, 0, 3, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 66, 0, 3, 7, 8,
                                                                       24, 27, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 72, 0, 3, 8, 9,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 78, 0, 3, 9, 10,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 84, 0, 3, 10, 11,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 90, 0, 3, 11, 12,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 96, 0, 3, 12, 13,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 102, 0, 3, 16, 17,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 108, 0, 3, 17, 18,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 114, 0, 3, 18, 19,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 120, 0, 3, 19, 20,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 126, 0, 3, 20, 21,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 132, 0, 3, 21, 22,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 138, 0, 3, 24, 27,
                                                                       66, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 148, 0, 3, 27, 30,
                                                                       72, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 158, 0, 3, 30, 33,
                                                                       78, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 168, 0, 3, 33, 36,
                                                                       84, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 178, 0, 3, 36, 39,
                                                                       90, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 188, 0, 3, 45, 48,
                                                                       102, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 198, 0, 3, 48, 51,
                                                                       108, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 208, 0, 3, 51, 54,
                                                                       114, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 218, 0, 3, 54, 57,
                                                                       120, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 228, 0, 3, 57, 60,
                                                                       126, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 238, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 241, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 244, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 247, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 250, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 253, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 256, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 259, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 262, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 265, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 268, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 271, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 274, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 277, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 280, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 283, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 286, 3, 7, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 295, 3, 8, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 304, 3, 9, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 313, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 322, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 331, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 340, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 349, 3, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 358, 3, 17, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 367, 3, 18, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 376, 3, 19, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 385, 3, 20, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 394, 3, 21, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 403, 3, 22, 63,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 412, 3, 24, 66,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 430, 3, 27, 72,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 448, 3, 30, 78,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 466, 3, 33, 84,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 484, 3, 36, 90,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 502, 3, 39, 96,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 520, 3, 45, 102,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 538, 3, 48, 108,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 556, 3, 51, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 574, 3, 54, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 592, 3, 57, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 610, 3, 60, 132,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 628, 3, 66, 138,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 658, 3, 72, 148,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 688, 3, 78, 158,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 718, 3, 84, 168,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 748, 3, 90, 178,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 778, 3, 102, 188,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 808, 3, 108, 198,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 838, 3, 114, 208,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 868, 3, 120, 218,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 898, 3, 126, 228,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 928, 3, 7, 8, 244,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 934, 3, 8, 9, 247,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 940, 3, 9, 10,
                                                                       250, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 946, 3, 10, 11,
                                                                       253, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 952, 3, 11, 12,
                                                                       256, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 958, 3, 12, 13,
                                                                       259, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 964, 3, 16, 17,
                                                                       268, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 970, 3, 17, 18,
                                                                       271, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 976, 3, 18, 19,
                                                                       274, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 982, 3, 19, 20,
                                                                       277, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 988, 3, 20, 21,
                                                                       280, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 994, 3, 21, 22,
                                                                       283, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1000, 0, 3, 928,
                                                                       244, 934, 304, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1018, 0, 3, 934,
                                                                       247, 940, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1036, 0, 3, 940,
                                                                       250, 946, 322, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1054, 0, 3, 946,
                                                                       253, 952, 331, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1072, 0, 3, 952,
                                                                       256, 958, 340, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1090, 0, 3, 964,
                                                                       268, 970, 367, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1108, 0, 3, 970,
                                                                       271, 976, 376, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1126, 0, 3, 976,
                                                                       274, 982, 385, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1144, 0, 3, 982,
                                                                       277, 988, 394, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1162, 0, 3, 988,
                                                                       280, 994, 403, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1180, 0, 3, 1000,
                                                                       304, 1018, 66, 72, 448,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1216, 0, 3, 1018,
                                                                       313, 1036, 72, 78, 466,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1252, 0, 3, 1036,
                                                                       322, 1054, 78, 84, 484,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1288, 0, 3, 1054,
                                                                       331, 1072, 84, 90, 502,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1324, 0, 3, 1090,
                                                                       367, 1108, 102, 108, 556,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1360, 0, 3, 1108,
                                                                       376, 1126, 108, 114, 574,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1396, 0, 3, 1126,
                                                                       385, 1144, 114, 120, 592,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1432, 0, 3, 1144,
                                                                       394, 1162, 120, 126, 610,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1468, 0, 3, 1180,
                                                                       448, 1216, 138, 148, 688,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1528, 0, 3, 1216,
                                                                       466, 1252, 148, 158, 718,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1588, 0, 3, 1252,
                                                                       484, 1288, 158, 168, 748,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1648, 0, 3, 1324,
                                                                       556, 1360, 188, 198, 838,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1708, 0, 3, 1360,
                                                                       574, 1396, 198, 208, 868,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1768, 0, 3, 1396,
                                                                       592, 1432, 208, 218, 898,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1828, 3, 238, 241,
                                                                       928, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1838, 3, 241, 244,
                                                                       934, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1848, 3, 244, 247,
                                                                       940, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1858, 3, 247, 250,
                                                                       946, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1868, 3, 250, 253,
                                                                       952, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1878, 3, 253, 256,
                                                                       958, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1888, 3, 262, 265,
                                                                       964, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1898, 3, 265, 268,
                                                                       970, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1908, 3, 268, 271,
                                                                       976, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1918, 3, 271, 274,
                                                                       982, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1928, 3, 274, 277,
                                                                       988, ncols, gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 1938, 3, 277, 280,
                                                                       994, ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1948, 0, 3, 1828,
                                                                       928, 1838, 286, 295, 1000,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 1978, 0, 3, 1838,
                                                                       934, 1848, 295, 304, 1018,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2008, 0, 3, 1848,
                                                                       940, 1858, 304, 313, 1036,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2038, 0, 3, 1858,
                                                                       946, 1868, 313, 322, 1054,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2068, 0, 3, 1868,
                                                                       952, 1878, 322, 331, 1072,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2098, 0, 3, 1888,
                                                                       964, 1898, 349, 358, 1090,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2128, 0, 3, 1898,
                                                                       970, 1908, 358, 367, 1108,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2158, 0, 3, 1908,
                                                                       976, 1918, 367, 376, 1126,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2188, 0, 3, 1918,
                                                                       982, 1928, 376, 385, 1144,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 2218, 0, 3, 1928,
                                                                       988, 1938, 385, 394, 1162,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2248, 0, 3, 1948,
                                                                       1000, 1978, 412, 430,
                                                                       1180, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2308, 0, 3, 1978,
                                                                       1018, 2008, 430, 448,
                                                                       1216, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2368, 0, 3, 2008,
                                                                       1036, 2038, 448, 466,
                                                                       1252, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2428, 0, 3, 2038,
                                                                       1054, 2068, 466, 484,
                                                                       1288, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2488, 0, 3, 2098,
                                                                       1090, 2128, 520, 538,
                                                                       1324, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2548, 0, 3, 2128,
                                                                       1108, 2158, 538, 556,
                                                                       1360, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2608, 0, 3, 2158,
                                                                       1126, 2188, 556, 574,
                                                                       1396, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 2668, 0, 3, 2188,
                                                                       1144, 2218, 574, 592,
                                                                       1432, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 2728, 0, 3, 2248,
                                                                       1180, 2308, 628, 658,
                                                                       1468, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 2828, 0, 3, 2308,
                                                                       1216, 2368, 658, 688,
                                                                       1528, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 2928, 0, 3, 2368,
                                                                       1252, 2428, 688, 718,
                                                                       1588, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3028, 0, 3, 2488,
                                                                       1324, 2548, 778, 808,
                                                                       1648, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3128, 0, 3, 2548,
                                                                       1360, 2608, 808, 838,
                                                                       1708, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 3228, 0, 3, 2608,
                                                                       1396, 2668, 838, 868,
                                                                       1768, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3328, 3, 928, 934,
                                                                       1848, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3343, 3, 934, 940,
                                                                       1858, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3358, 3, 940, 946,
                                                                       1868, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3373, 3, 946, 952,
                                                                       1878, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3388, 3, 964, 970,
                                                                       1908, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3403, 3, 970, 976,
                                                                       1918, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3418, 3, 976, 982,
                                                                       1928, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 3433, 3, 982, 988,
                                                                       1938, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3448, 0, 3, 3328,
                                                                       1848, 3343, 1000, 1018,
                                                                       2008, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3493, 0, 3, 3343,
                                                                       1858, 3358, 1018, 1036,
                                                                       2038, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3538, 0, 3, 3358,
                                                                       1868, 3373, 1036, 1054,
                                                                       2068, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3583, 0, 3, 3388,
                                                                       1908, 3403, 1090, 1108,
                                                                       2158, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3628, 0, 3, 3403,
                                                                       1918, 3418, 1108, 1126,
                                                                       2188, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 3673, 0, 3, 3418,
                                                                       1928, 3433, 1126, 1144,
                                                                       2218, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 3718, 0, 3, 3448,
                                                                       2008, 3493, 1180, 1216,
                                                                       2368, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 3808, 0, 3, 3493,
                                                                       2038, 3538, 1216, 1252,
                                                                       2428, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 3898, 0, 3, 3583,
                                                                       2158, 3628, 1324, 1360,
                                                                       2608, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 3988, 0, 3, 3628,
                                                                       2188, 3673, 1360, 1396,
                                                                       2668, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 4078, 0, 3, 3718,
                                                                       2368, 3808, 1468, 1528,
                                                                       2928, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 4228, 0, 3, 3898,
                                                                       2608, 3988, 1648, 1708,
                                                                       3228, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4378, 3, 1828,
                                                                       1838, 3328, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4399, 3, 1838,
                                                                       1848, 3343, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4420, 3, 1848,
                                                                       1858, 3358, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4441, 3, 1858,
                                                                       1868, 3373, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4462, 3, 1888,
                                                                       1898, 3388, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4483, 3, 1898,
                                                                       1908, 3403, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4504, 3, 1908,
                                                                       1918, 3418, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 4525, 3, 1918,
                                                                       1928, 3433, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 4546, 0, 3, 4378,
                                                                       3328, 4399, 1948, 1978,
                                                                       3448, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 4609, 0, 3, 4399,
                                                                       3343, 4420, 1978, 2008,
                                                                       3493, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 4672, 0, 3, 4420,
                                                                       3358, 4441, 2008, 2038,
                                                                       3538, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 4735, 0, 3, 4462,
                                                                       3388, 4483, 2098, 2128,
                                                                       3583, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 4798, 0, 3, 4483,
                                                                       3403, 4504, 2128, 2158,
                                                                       3628, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 4861, 0, 3, 4504,
                                                                       3418, 4525, 2158, 2188,
                                                                       3673, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 4924, 0, 3, 4546,
                                                                       3448, 4609, 2248, 2308,
                                                                       3718, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 5050, 0, 3, 4609,
                                                                       3493, 4672, 2308, 2368,
                                                                       3808, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 5176, 0, 3, 4735,
                                                                       3583, 4798, 2488, 2548,
                                                                       3898, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 5302, 0, 3, 4798,
                                                                       3628, 4861, 2548, 2608,
                                                                       3988, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 5428, 0, 3, 4924,
                                                                       3718, 5050, 2728, 2828,
                                                                       4078, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 5638, 0, 3, 5176,
                                                                       3898, 5302, 3028, 3128,
                                                                       4228, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_p_x(buffer, 5848, 4378, 4924, 1, 21, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 5911, 4378, 4924, 1, 21, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 5974, 4378, 4924, 1, 21, ncols, alpha);

                    simdgeo::geom_p_x(buffer, 6037, 4462, 5176, 1, 21, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 6100, 4462, 5176, 1, 21, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 6163, 4462, 5176, 1, 21, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 6226, 4546, 5428, 1, 21, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 6352, 4546, 5428, 1, 21, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 6478, 4546, 5428, 1, 21, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 6604, 4735, 5638, 1, 21, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 6730, 4735, 5638, 1, 21, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 6856, 4735, 5638, 1, 21, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 6982, 5848, 63, ncols);

                    simdfunc::contract_primitives(buffer, 7078, 5911, 63, ncols);

                    simdfunc::contract_primitives(buffer, 7174, 5974, 63, ncols);

                    simdfunc::contract_primitives(buffer, 7270, 4546, 63, ncols);

                    simdfunc::contract_primitives(buffer, 7366, 6037, 63, ncols);

                    simdfunc::contract_primitives(buffer, 7462, 6100, 63, ncols);

                    simdfunc::contract_primitives(buffer, 7558, 6163, 63, ncols);

                    simdfunc::contract_primitives(buffer, 7654, 4735, 63, ncols);

                    simdfunc::contract_primitives(buffer, 7750, 6226, 126, ncols);

                    simdfunc::contract_primitives(buffer, 7942, 6352, 126, ncols);

                    simdfunc::contract_primitives(buffer, 8134, 6478, 126, ncols);

                    simdfunc::contract_primitives(buffer, 8326, 6604, 126, ncols);

                    simdfunc::contract_primitives(buffer, 8518, 6730, 126, ncols);

                    simdfunc::contract_primitives(buffer, 8710, 6856, 126, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 7045, 6982, 3, 1, nmax);

        simdtrf::transform_h_inner(buffer, 7141, 7078, 3, 1, nmax);

        simdtrf::transform_h_inner(buffer, 7237, 7174, 3, 1, nmax);

        simdtrf::transform_h_inner(buffer, 7333, 7270, 3, 1, nmax);

        simdtrf::transform_h_inner(buffer, 7429, 7366, 3, 1, nmax);

        simdtrf::transform_h_inner(buffer, 7525, 7462, 3, 1, nmax);

        simdtrf::transform_h_inner(buffer, 7621, 7558, 3, 1, nmax);

        simdtrf::transform_h_inner(buffer, 7717, 7654, 3, 1, nmax);

        simdtrf::transform_h_inner(buffer, 7876, 7750, 6, 1, nmax);

        simdtrf::transform_h_inner(buffer, 8068, 7942, 6, 1, nmax);

        simdtrf::transform_h_inner(buffer, 8260, 8134, 6, 1, nmax);

        simdtrf::transform_h_inner(buffer, 8452, 8326, 6, 1, nmax);

        simdtrf::transform_h_inner(buffer, 8644, 8518, 6, 1, nmax);

        simdtrf::transform_h_inner(buffer, 8836, 8710, 6, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 8902, 7045, 7333,
                                                       7876, 11, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 9001, 7141, 7333,
                                                       8068, 11, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 9100, 7237, 7333,
                                                       8260, 11, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 9199, 7429, 7717,
                                                       8452, 11, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 9298, 7525, 7717,
                                                       8644, 11, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 9397, 7621, 7717,
                                                       8836, 11, nmax);

        simdtrf::transform_p_inner(buffer, 9496, 9199, 3, 11, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 9496, 33, nmax);

        simdtrf::transform_p_inner(buffer, 9496, 9298, 3, 11, nmax);

        simdtrf::transform_p_outer(values + 99 * nvalues + n * npairs, nvalues, buffer, 9496, 33,
                                   nmax);

        simdtrf::transform_p_inner(buffer, 9496, 9397, 3, 11, nmax);

        simdtrf::transform_p_outer(values + 198 * nvalues + n * npairs, nvalues, buffer, 9496,
                                   33, nmax);

        simdtrf::transform_p_inner(buffer, 9496, 8902, 3, 11, nmax);

        simdtrf::transform_p_outer(values + 297 * nvalues + n * npairs, nvalues, buffer, 9496,
                                   33, nmax);

        simdtrf::transform_p_inner(buffer, 9496, 9001, 3, 11, nmax);

        simdtrf::transform_p_outer(values + 396 * nvalues + n * npairs, nvalues, buffer, 9496,
                                   33, nmax);

        simdtrf::transform_p_inner(buffer, 9496, 9100, 3, 11, nmax);

        simdtrf::transform_p_outer(values + 495 * nvalues + n * npairs, nvalues, buffer, 9496,
                                   33, nmax);
    }

    for (size_t m = 0; m < 594; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
