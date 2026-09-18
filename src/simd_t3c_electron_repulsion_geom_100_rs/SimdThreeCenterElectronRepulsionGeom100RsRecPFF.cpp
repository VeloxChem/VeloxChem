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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecPFF.hpp"

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
#include "SimdGeometryP1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecHSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferGeom100XDD.hpp"
#include "SimdTransferGeom100XDP.hpp"
#include "SimdTransferGeom100XFP.hpp"
#include "SimdTransferGeom100XPD.hpp"
#include "SimdTransferGeom100XPF.hpp"
#include "SimdTransferGeom100XPP.hpp"
#include "SimdTransferGeom100YDD.hpp"
#include "SimdTransferGeom100YDP.hpp"
#include "SimdTransferGeom100YFP.hpp"
#include "SimdTransferGeom100YPD.hpp"
#include "SimdTransferGeom100YPF.hpp"
#include "SimdTransferGeom100YPP.hpp"
#include "SimdTransferGeom100ZDD.hpp"
#include "SimdTransferGeom100ZDP.hpp"
#include "SimdTransferGeom100ZFP.hpp"
#include "SimdTransferGeom100ZPD.hpp"
#include "SimdTransferGeom100ZPF.hpp"
#include "SimdTransferGeom100ZPP.hpp"
#include "SimdTransferPD.hpp"
#include "SimdTransferPP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformP.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_pff_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_pff_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 18797, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 882 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 18797, 7984, 4009, dimensions);

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

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 238, 0, 3, 66, 72,
                                                                       138, 148, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 253, 0, 3, 72, 78,
                                                                       148, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 268, 0, 3, 78, 84,
                                                                       158, 168, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 283, 0, 3, 84, 90,
                                                                       168, 178, ncols, gamma, p,
                                                                       q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 298, 0, 3, 102,
                                                                       108, 188, 198, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 313, 0, 3, 108,
                                                                       114, 198, 208, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 328, 0, 3, 114,
                                                                       120, 208, 218, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 343, 0, 3, 120,
                                                                       126, 218, 228, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 358, 0, 3, 138,
                                                                       148, 238, 253, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 379, 0, 3, 148,
                                                                       158, 253, 268, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 400, 0, 3, 158,
                                                                       168, 268, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 421, 0, 3, 188,
                                                                       198, 298, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 442, 0, 3, 198,
                                                                       208, 313, 328, ncols,
                                                                       gamma, p, q);

                    compute_prim_hss_three_center_electron_repulsion_0(buffer, 463, 0, 3, 208,
                                                                       218, 328, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 484, 3, 7, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 487, 3, 8, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 490, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 493, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 496, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 499, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 502, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 505, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 508, 3, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 511, 3, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 514, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 517, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 520, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 523, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 526, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 529, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 532, 3, 7, 24,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 541, 3, 8, 27,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 550, 3, 9, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 559, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 568, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 577, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 586, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 595, 3, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 604, 3, 17, 48,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 613, 3, 18, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 622, 3, 19, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 631, 3, 20, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 640, 3, 21, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 649, 3, 22, 63,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 658, 3, 24, 66,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 676, 3, 27, 72,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 694, 3, 30, 78,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 712, 3, 33, 84,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 730, 3, 36, 90,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 748, 3, 39, 96,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 766, 3, 45, 102,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 784, 3, 48, 108,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 802, 3, 51, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 820, 3, 54, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 838, 3, 57, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 856, 3, 60, 132,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 874, 3, 66, 138,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 904, 3, 72, 148,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 934, 3, 78, 158,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 964, 3, 84, 168,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 994, 3, 90, 178,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1024, 3, 102, 188,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1054, 3, 108, 198,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1084, 3, 114, 208,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1114, 3, 120, 218,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1144, 3, 126, 228,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1174, 3, 138, 238,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1219, 3, 148, 253,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1264, 3, 158, 268,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1309, 3, 168, 283,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1354, 3, 188, 298,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1399, 3, 198, 313,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1444, 3, 208, 328,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1489, 3, 218, 343,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1534, 3, 238, 358,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1597, 3, 253, 379,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1660, 3, 268, 400,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1723, 3, 298, 421,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1786, 3, 313, 442,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1849, 3, 328, 463,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1912, 3, 7, 8,
                                                                       490, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1918, 3, 8, 9,
                                                                       493, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1924, 3, 9, 10,
                                                                       496, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1930, 3, 10, 11,
                                                                       499, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1936, 3, 11, 12,
                                                                       502, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1942, 3, 12, 13,
                                                                       505, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1948, 3, 16, 17,
                                                                       514, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1954, 3, 17, 18,
                                                                       517, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1960, 3, 18, 19,
                                                                       520, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1966, 3, 19, 20,
                                                                       523, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1972, 3, 20, 21,
                                                                       526, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1978, 3, 21, 22,
                                                                       529, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1984, 0, 3, 1912,
                                                                       490, 1918, 550, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2002, 0, 3, 1918,
                                                                       493, 1924, 559, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2020, 0, 3, 1924,
                                                                       496, 1930, 568, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2038, 0, 3, 1930,
                                                                       499, 1936, 577, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2056, 0, 3, 1936,
                                                                       502, 1942, 586, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2074, 0, 3, 1948,
                                                                       514, 1954, 613, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2092, 0, 3, 1954,
                                                                       517, 1960, 622, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2110, 0, 3, 1960,
                                                                       520, 1966, 631, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2128, 0, 3, 1966,
                                                                       523, 1972, 640, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2146, 0, 3, 1972,
                                                                       526, 1978, 649, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2164, 0, 3, 1984,
                                                                       550, 2002, 66, 72, 694,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2200, 0, 3, 2002,
                                                                       559, 2020, 72, 78, 712,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2236, 0, 3, 2020,
                                                                       568, 2038, 78, 84, 730,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2272, 0, 3, 2038,
                                                                       577, 2056, 84, 90, 748,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2308, 0, 3, 2074,
                                                                       613, 2092, 102, 108, 802,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2344, 0, 3, 2092,
                                                                       622, 2110, 108, 114, 820,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2380, 0, 3, 2110,
                                                                       631, 2128, 114, 120, 838,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 2416, 0, 3, 2128,
                                                                       640, 2146, 120, 126, 856,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2452, 0, 3, 2164,
                                                                       694, 2200, 138, 148, 934,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2512, 0, 3, 2200,
                                                                       712, 2236, 148, 158, 964,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2572, 0, 3, 2236,
                                                                       730, 2272, 158, 168, 994,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2632, 0, 3, 2308,
                                                                       802, 2344, 188, 198, 1084,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2692, 0, 3, 2344,
                                                                       820, 2380, 198, 208, 1114,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2752, 0, 3, 2380,
                                                                       838, 2416, 208, 218, 1144,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2812, 0, 3, 2452,
                                                                       934, 2512, 238, 253, 1264,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2902, 0, 3, 2512,
                                                                       964, 2572, 253, 268, 1309,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2992, 0, 3, 2632,
                                                                       1084, 2692, 298, 313,
                                                                       1444, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 3082, 0, 3, 2692,
                                                                       1114, 2752, 313, 328,
                                                                       1489, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 3172, 0, 3, 2812,
                                                                       1264, 2902, 358, 379,
                                                                       1660, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 3298, 0, 3, 2992,
                                                                       1444, 3082, 421, 442,
                                                                       1849, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3424, 3, 484, 487,
                                                                       1912, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3434, 3, 487, 490,
                                                                       1918, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3444, 3, 490, 493,
                                                                       1924, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3454, 3, 493, 496,
                                                                       1930, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3464, 3, 496, 499,
                                                                       1936, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3474, 3, 499, 502,
                                                                       1942, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3484, 3, 508, 511,
                                                                       1948, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3494, 3, 511, 514,
                                                                       1954, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3504, 3, 514, 517,
                                                                       1960, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3514, 3, 517, 520,
                                                                       1966, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3524, 3, 520, 523,
                                                                       1972, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3534, 3, 523, 526,
                                                                       1978, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3544, 0, 3, 3424,
                                                                       1912, 3434, 532, 541,
                                                                       1984, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3574, 0, 3, 3434,
                                                                       1918, 3444, 541, 550,
                                                                       2002, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3604, 0, 3, 3444,
                                                                       1924, 3454, 550, 559,
                                                                       2020, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3634, 0, 3, 3454,
                                                                       1930, 3464, 559, 568,
                                                                       2038, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3664, 0, 3, 3464,
                                                                       1936, 3474, 568, 577,
                                                                       2056, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3694, 0, 3, 3484,
                                                                       1948, 3494, 595, 604,
                                                                       2074, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3724, 0, 3, 3494,
                                                                       1954, 3504, 604, 613,
                                                                       2092, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3754, 0, 3, 3504,
                                                                       1960, 3514, 613, 622,
                                                                       2110, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3784, 0, 3, 3514,
                                                                       1966, 3524, 622, 631,
                                                                       2128, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 3814, 0, 3, 3524,
                                                                       1972, 3534, 631, 640,
                                                                       2146, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3844, 0, 3, 3544,
                                                                       1984, 3574, 658, 676,
                                                                       2164, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3904, 0, 3, 3574,
                                                                       2002, 3604, 676, 694,
                                                                       2200, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 3964, 0, 3, 3604,
                                                                       2020, 3634, 694, 712,
                                                                       2236, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4024, 0, 3, 3634,
                                                                       2038, 3664, 712, 730,
                                                                       2272, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4084, 0, 3, 3694,
                                                                       2074, 3724, 766, 784,
                                                                       2308, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4144, 0, 3, 3724,
                                                                       2092, 3754, 784, 802,
                                                                       2344, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4204, 0, 3, 3754,
                                                                       2110, 3784, 802, 820,
                                                                       2380, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 4264, 0, 3, 3784,
                                                                       2128, 3814, 820, 838,
                                                                       2416, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4324, 0, 3, 3844,
                                                                       2164, 3904, 874, 904,
                                                                       2452, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4424, 0, 3, 3904,
                                                                       2200, 3964, 904, 934,
                                                                       2512, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4524, 0, 3, 3964,
                                                                       2236, 4024, 934, 964,
                                                                       2572, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4624, 0, 3, 4084,
                                                                       2308, 4144, 1024, 1054,
                                                                       2632, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4724, 0, 3, 4144,
                                                                       2344, 4204, 1054, 1084,
                                                                       2692, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 4824, 0, 3, 4204,
                                                                       2380, 4264, 1084, 1114,
                                                                       2752, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 4924, 0, 3, 4324,
                                                                       2452, 4424, 1174, 1219,
                                                                       2812, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 5074, 0, 3, 4424,
                                                                       2512, 4524, 1219, 1264,
                                                                       2902, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 5224, 0, 3, 4624,
                                                                       2632, 4724, 1354, 1399,
                                                                       2992, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 5374, 0, 3, 4724,
                                                                       2692, 4824, 1399, 1444,
                                                                       3082, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 5524, 0, 3, 4924,
                                                                       2812, 5074, 1534, 1597,
                                                                       3172, ncols, gamma, p,
                                                                       q);

                    compute_prim_hsf_three_center_electron_repulsion_0(buffer, 5734, 0, 3, 5224,
                                                                       2992, 5374, 1723, 1786,
                                                                       3298, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_p_x(buffer, 5944, 3424, 3844, 1, 10, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 5974, 3424, 3844, 1, 10, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 6004, 3424, 3844, 1, 10, ncols, alpha);

                    simdgeo::geom_p_x(buffer, 6034, 3484, 4084, 1, 10, ncols, alpha);

                    simdgeo::geom_p_y(buffer, 6064, 3484, 4084, 1, 10, ncols, alpha);

                    simdgeo::geom_p_z(buffer, 6094, 3484, 4084, 1, 10, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 6124, 3544, 4324, 1, 10, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 6184, 3544, 4324, 1, 10, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 6244, 3544, 4324, 1, 10, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 6304, 3694, 4624, 1, 10, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 6364, 3694, 4624, 1, 10, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 6424, 3694, 4624, 1, 10, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 6484, 3844, 4924, 1, 10, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 6584, 3844, 4924, 1, 10, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 6684, 3844, 4924, 1, 10, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 6784, 4084, 5224, 1, 10, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 6884, 4084, 5224, 1, 10, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 6984, 4084, 5224, 1, 10, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 7084, 4324, 5524, 1, 10, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 7234, 4324, 5524, 1, 10, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 7384, 4324, 5524, 1, 10, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 7534, 4624, 5734, 1, 10, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 7684, 4624, 5734, 1, 10, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 7834, 4624, 5734, 1, 10, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 7984, 5944, 30, ncols);

                    simdfunc::contract_primitives(buffer, 8035, 5974, 30, ncols);

                    simdfunc::contract_primitives(buffer, 8086, 6004, 30, ncols);

                    simdfunc::contract_primitives(buffer, 8137, 3544, 30, ncols);

                    simdfunc::contract_primitives(buffer, 8188, 6034, 30, ncols);

                    simdfunc::contract_primitives(buffer, 8239, 6064, 30, ncols);

                    simdfunc::contract_primitives(buffer, 8290, 6094, 30, ncols);

                    simdfunc::contract_primitives(buffer, 8341, 3694, 30, ncols);

                    simdfunc::contract_primitives(buffer, 8392, 6124, 60, ncols);

                    simdfunc::contract_primitives(buffer, 8494, 6184, 60, ncols);

                    simdfunc::contract_primitives(buffer, 8596, 6244, 60, ncols);

                    simdfunc::contract_primitives(buffer, 8698, 3844, 60, ncols);

                    simdfunc::contract_primitives(buffer, 8800, 6304, 60, ncols);

                    simdfunc::contract_primitives(buffer, 8902, 6364, 60, ncols);

                    simdfunc::contract_primitives(buffer, 9004, 6424, 60, ncols);

                    simdfunc::contract_primitives(buffer, 9106, 4084, 60, ncols);

                    simdfunc::contract_primitives(buffer, 9208, 6484, 100, ncols);

                    simdfunc::contract_primitives(buffer, 9378, 6584, 100, ncols);

                    simdfunc::contract_primitives(buffer, 9548, 6684, 100, ncols);

                    simdfunc::contract_primitives(buffer, 9718, 4324, 100, ncols);

                    simdfunc::contract_primitives(buffer, 9888, 6784, 100, ncols);

                    simdfunc::contract_primitives(buffer, 10058, 6884, 100, ncols);

                    simdfunc::contract_primitives(buffer, 10228, 6984, 100, ncols);

                    simdfunc::contract_primitives(buffer, 10398, 4624, 100, ncols);

                    simdfunc::contract_primitives(buffer, 10568, 7084, 150, ncols);

                    simdfunc::contract_primitives(buffer, 10823, 7234, 150, ncols);

                    simdfunc::contract_primitives(buffer, 11078, 7384, 150, ncols);

                    simdfunc::contract_primitives(buffer, 11333, 7534, 150, ncols);

                    simdfunc::contract_primitives(buffer, 11588, 7684, 150, ncols);

                    simdfunc::contract_primitives(buffer, 11843, 7834, 150, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 8014, 7984, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8065, 8035, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8116, 8086, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8167, 8137, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8218, 8188, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8269, 8239, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8320, 8290, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8371, 8341, 3, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8452, 8392, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8554, 8494, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8656, 8596, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8758, 8698, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8860, 8800, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 8962, 8902, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 9064, 9004, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 9166, 9106, 6, 1, nmax);

        simdtrf::transform_f_inner(buffer, 9308, 9208, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 9478, 9378, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 9648, 9548, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 9818, 9718, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 9988, 9888, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10158, 10058, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10328, 10228, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10498, 10398, 10, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10718, 10568, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 10973, 10823, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 11228, 11078, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 11483, 11333, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 11738, 11588, 15, 1, nmax);

        simdtrf::transform_f_inner(buffer, 11993, 11843, 15, 1, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 12098, 8014, 8167,
                                                       8452, 7, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 12161, 8065, 8167,
                                                       8554, 7, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 12224, 8116, 8167,
                                                       8656, 7, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 12287, 8167, 8758, 7, nmax);

        simdtrf::compute_hrr_geom_100x_pp_out_of_first(buffer, coordinates, 12350, 8218, 8371,
                                                       8860, 7, nmax);

        simdtrf::compute_hrr_geom_100y_pp_out_of_first(buffer, coordinates, 12413, 8269, 8371,
                                                       8962, 7, nmax);

        simdtrf::compute_hrr_geom_100z_pp_out_of_first(buffer, coordinates, 12476, 8320, 8371,
                                                       9064, 7, nmax);

        simdtrf::compute_hrr_pp_out_of_first(buffer, coordinates, 12539, 8371, 9166, 7, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 12602, 8452, 8758,
                                                       9308, 7, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 12728, 8554, 8758,
                                                       9478, 7, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 12854, 8656, 8758,
                                                       9648, 7, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 12980, 8758, 9818, 7, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 13106, 8860, 9166,
                                                       9988, 7, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 13232, 8962, 9166,
                                                       10158, 7, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 13358, 9064, 9166,
                                                       10328, 7, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 13484, 9166, 10498, 7, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 13610, 9308, 9818,
                                                       10718, 7, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 13820, 9478, 9818,
                                                       10973, 7, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 14030, 9648, 9818,
                                                       11228, 7, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 14240, 9988, 10498,
                                                       11483, 7, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 14450, 10158, 10498,
                                                       11738, 7, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 14660, 10328, 10498,
                                                       11993, 7, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 14870, 12098, 12287,
                                                       12602, 7, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 14996, 12161, 12287,
                                                       12728, 7, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 15122, 12224, 12287,
                                                       12854, 7, nmax);

        simdtrf::compute_hrr_pd_out_of_first(buffer, coordinates, 15248, 12287, 12980, 7, nmax);

        simdtrf::compute_hrr_geom_100x_pd_out_of_first(buffer, coordinates, 15374, 12350, 12539,
                                                       13106, 7, nmax);

        simdtrf::compute_hrr_geom_100y_pd_out_of_first(buffer, coordinates, 15500, 12413, 12539,
                                                       13232, 7, nmax);

        simdtrf::compute_hrr_geom_100z_pd_out_of_first(buffer, coordinates, 15626, 12476, 12539,
                                                       13358, 7, nmax);

        simdtrf::compute_hrr_pd_out_of_first(buffer, coordinates, 15752, 12539, 13484, 7, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 15878, 12602, 12980,
                                                       13610, 7, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 16130, 12728, 12980,
                                                       13820, 7, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 16382, 12854, 12980,
                                                       14030, 7, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 16634, 13106, 13484,
                                                       14240, 7, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 16886, 13232, 13484,
                                                       14450, 7, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 17138, 13358, 13484,
                                                       14660, 7, nmax);

        simdtrf::compute_hrr_geom_100x_pf_out_of_first(buffer, coordinates, 17390, 14870, 15248,
                                                       15878, 7, nmax);

        simdtrf::compute_hrr_geom_100y_pf_out_of_first(buffer, coordinates, 17600, 14996, 15248,
                                                       16130, 7, nmax);

        simdtrf::compute_hrr_geom_100z_pf_out_of_first(buffer, coordinates, 17810, 15122, 15248,
                                                       16382, 7, nmax);

        simdtrf::compute_hrr_geom_100x_pf_out_of_first(buffer, coordinates, 18020, 15374, 15752,
                                                       16634, 7, nmax);

        simdtrf::compute_hrr_geom_100y_pf_out_of_first(buffer, coordinates, 18230, 15500, 15752,
                                                       16886, 7, nmax);

        simdtrf::compute_hrr_geom_100z_pf_out_of_first(buffer, coordinates, 18440, 15626, 15752,
                                                       17138, 7, nmax);

        simdtrf::transform_f_inner(buffer, 18650, 18020, 3, 7, nmax);

        simdtrf::transform_p_outer(values + n * npairs, nvalues, buffer, 18650, 49, nmax);

        simdtrf::transform_f_inner(buffer, 18650, 18230, 3, 7, nmax);

        simdtrf::transform_p_outer(values + 147 * nvalues + n * npairs, nvalues, buffer, 18650,
                                   49, nmax);

        simdtrf::transform_f_inner(buffer, 18650, 18440, 3, 7, nmax);

        simdtrf::transform_p_outer(values + 294 * nvalues + n * npairs, nvalues, buffer, 18650,
                                   49, nmax);

        simdtrf::transform_f_inner(buffer, 18650, 17390, 3, 7, nmax);

        simdtrf::transform_p_outer(values + 441 * nvalues + n * npairs, nvalues, buffer, 18650,
                                   49, nmax);

        simdtrf::transform_f_inner(buffer, 18650, 17600, 3, 7, nmax);

        simdtrf::transform_p_outer(values + 588 * nvalues + n * npairs, nvalues, buffer, 18650,
                                   49, nmax);

        simdtrf::transform_f_inner(buffer, 18650, 17810, 3, 7, nmax);

        simdtrf::transform_p_outer(values + 735 * nvalues + n * npairs, nvalues, buffer, 18650,
                                   49, nmax);
    }

    for (size_t m = 0; m < 882; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
