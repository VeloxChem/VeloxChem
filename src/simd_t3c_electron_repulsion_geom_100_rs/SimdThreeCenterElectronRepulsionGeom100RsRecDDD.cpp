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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecDDD.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransferDP.hpp"
#include "SimdTransferGeom100XDD.hpp"
#include "SimdTransferGeom100XDP.hpp"
#include "SimdTransferGeom100XFP.hpp"
#include "SimdTransferGeom100YDD.hpp"
#include "SimdTransferGeom100YDP.hpp"
#include "SimdTransferGeom100YFP.hpp"
#include "SimdTransferGeom100ZDD.hpp"
#include "SimdTransferGeom100ZDP.hpp"
#include "SimdTransferGeom100ZFP.hpp"
#include "SimdTransformD.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_ddd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_ddd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 9116, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 750 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 9116, 3868, 2323, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 6, 3, 7,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 15, 3, 7,
                                                             ncols, fj, i * nprim_b + j, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 484, 3, 9, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 487, 3, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 490, 3, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 493, 3, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 496, 3, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 499, 3, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 502, 3, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 505, 3, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 508, 3, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 511, 3, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 514, 3, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 517, 3, 23, ncols,
                                                                       p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 520, 3, 9, 30,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 529, 3, 10, 33,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 538, 3, 11, 36,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 547, 3, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 556, 3, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 565, 3, 18, 51,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 574, 3, 19, 54,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 583, 3, 20, 57,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 592, 3, 21, 60,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 601, 3, 22, 63,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 610, 3, 30, 78,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 628, 3, 33, 84,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 646, 3, 36, 90,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 664, 3, 39, 96,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 682, 3, 51, 114,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 700, 3, 54, 120,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 718, 3, 57, 126,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 736, 3, 60, 132,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 754, 3, 78, 158,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 784, 3, 84, 168,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 814, 3, 90, 178,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 844, 3, 114, 208,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 874, 3, 120, 218,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 904, 3, 126, 228,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 934, 3, 158, 268,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 979, 3, 168, 283,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1024, 3, 208, 328,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 1069, 3, 218, 343,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1114, 3, 268, 400,
                                                                       ncols, p, q);

                    compute_prim_hsp_three_center_electron_repulsion_0(buffer, 1177, 3, 328, 463,
                                                                       ncols, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1240, 3, 7, 8,
                                                                       484, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1246, 3, 8, 9,
                                                                       487, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1252, 3, 9, 10,
                                                                       490, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1258, 3, 10, 11,
                                                                       493, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1264, 3, 11, 12,
                                                                       496, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1270, 3, 12, 13,
                                                                       499, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1276, 3, 16, 17,
                                                                       502, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1282, 3, 17, 18,
                                                                       505, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1288, 3, 18, 19,
                                                                       508, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1294, 3, 19, 20,
                                                                       511, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1300, 3, 20, 21,
                                                                       514, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1306, 3, 21, 22,
                                                                       517, ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1312, 0, 3, 1240,
                                                                       484, 1246, 520, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1330, 0, 3, 1246,
                                                                       487, 1252, 529, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1348, 0, 3, 1252,
                                                                       490, 1258, 538, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1366, 0, 3, 1258,
                                                                       493, 1264, 547, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1384, 0, 3, 1264,
                                                                       496, 1270, 556, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1402, 0, 3, 1276,
                                                                       502, 1282, 565, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1420, 0, 3, 1282,
                                                                       505, 1288, 574, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1438, 0, 3, 1288,
                                                                       508, 1294, 583, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1456, 0, 3, 1294,
                                                                       511, 1300, 592, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 1474, 0, 3, 1300,
                                                                       514, 1306, 601, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1492, 0, 3, 1312,
                                                                       520, 1330, 66, 72, 610,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1528, 0, 3, 1330,
                                                                       529, 1348, 72, 78, 628,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1564, 0, 3, 1348,
                                                                       538, 1366, 78, 84, 646,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1600, 0, 3, 1366,
                                                                       547, 1384, 84, 90, 664,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1636, 0, 3, 1402,
                                                                       565, 1420, 102, 108, 682,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1672, 0, 3, 1420,
                                                                       574, 1438, 108, 114, 700,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1708, 0, 3, 1438,
                                                                       583, 1456, 114, 120, 718,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1744, 0, 3, 1456,
                                                                       592, 1474, 120, 126, 736,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1780, 0, 3, 1492,
                                                                       610, 1528, 138, 148, 754,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1840, 0, 3, 1528,
                                                                       628, 1564, 148, 158, 784,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1900, 0, 3, 1564,
                                                                       646, 1600, 158, 168, 814,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 1960, 0, 3, 1636,
                                                                       682, 1672, 188, 198, 844,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2020, 0, 3, 1672,
                                                                       700, 1708, 198, 208, 874,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 2080, 0, 3, 1708,
                                                                       718, 1744, 208, 218, 904,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2140, 0, 3, 1780,
                                                                       754, 1840, 238, 253, 934,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2230, 0, 3, 1840,
                                                                       784, 1900, 253, 268, 979,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2320, 0, 3, 1960,
                                                                       844, 2020, 298, 313, 1024,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 2410, 0, 3, 2020,
                                                                       874, 2080, 313, 328, 1069,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 2500, 0, 3, 2140,
                                                                       934, 2230, 358, 379, 1114,
                                                                       ncols, gamma, p, q);

                    compute_prim_hsd_three_center_electron_repulsion_0(buffer, 2626, 0, 3, 2320,
                                                                       1024, 2410, 421, 442,
                                                                       1177, ncols, gamma, p,
                                                                       q);

                    simdgeo::geom_d_x(buffer, 2752, 1312, 1780, 1, 6, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 2788, 1312, 1780, 1, 6, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 2824, 1312, 1780, 1, 6, ncols, alpha);

                    simdgeo::geom_d_x(buffer, 2860, 1402, 1960, 1, 6, ncols, alpha);

                    simdgeo::geom_d_y(buffer, 2896, 1402, 1960, 1, 6, ncols, alpha);

                    simdgeo::geom_d_z(buffer, 2932, 1402, 1960, 1, 6, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 2968, 1492, 2140, 1, 6, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 3028, 1492, 2140, 1, 6, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 3088, 1492, 2140, 1, 6, ncols, alpha);

                    simdgeo::geom_f_x(buffer, 3148, 1636, 2320, 1, 6, ncols, alpha);

                    simdgeo::geom_f_y(buffer, 3208, 1636, 2320, 1, 6, ncols, alpha);

                    simdgeo::geom_f_z(buffer, 3268, 1636, 2320, 1, 6, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 3328, 1780, 2500, 1, 6, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 3418, 1780, 2500, 1, 6, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 3508, 1780, 2500, 1, 6, ncols, alpha);

                    simdgeo::geom_g_x(buffer, 3598, 1960, 2626, 1, 6, ncols, alpha);

                    simdgeo::geom_g_y(buffer, 3688, 1960, 2626, 1, 6, ncols, alpha);

                    simdgeo::geom_g_z(buffer, 3778, 1960, 2626, 1, 6, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 3868, 2752, 36, ncols);

                    simdfunc::contract_primitives(buffer, 3934, 2788, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4000, 2824, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4066, 1492, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4132, 2860, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4198, 2896, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4264, 2932, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4330, 1636, 36, ncols);

                    simdfunc::contract_primitives(buffer, 4396, 2968, 60, ncols);

                    simdfunc::contract_primitives(buffer, 4506, 3028, 60, ncols);

                    simdfunc::contract_primitives(buffer, 4616, 3088, 60, ncols);

                    simdfunc::contract_primitives(buffer, 4726, 1780, 60, ncols);

                    simdfunc::contract_primitives(buffer, 4836, 3148, 60, ncols);

                    simdfunc::contract_primitives(buffer, 4946, 3208, 60, ncols);

                    simdfunc::contract_primitives(buffer, 5056, 3268, 60, ncols);

                    simdfunc::contract_primitives(buffer, 5166, 1960, 60, ncols);

                    simdfunc::contract_primitives(buffer, 5276, 3328, 90, ncols);

                    simdfunc::contract_primitives(buffer, 5441, 3418, 90, ncols);

                    simdfunc::contract_primitives(buffer, 5606, 3508, 90, ncols);

                    simdfunc::contract_primitives(buffer, 5771, 3598, 90, ncols);

                    simdfunc::contract_primitives(buffer, 5936, 3688, 90, ncols);

                    simdfunc::contract_primitives(buffer, 6101, 3778, 90, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 3904, 3868, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 3970, 3934, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4036, 4000, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4102, 4066, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4168, 4132, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4234, 4198, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4300, 4264, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4366, 4330, 6, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4456, 4396, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4566, 4506, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4676, 4616, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4786, 4726, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 4896, 4836, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5006, 4946, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5116, 5056, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5226, 5166, 10, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5366, 5276, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5531, 5441, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5696, 5606, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 5861, 5771, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 6026, 5936, 15, 1, nmax);

        simdtrf::transform_d_inner(buffer, 6191, 6101, 15, 1, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 6266, 3904, 4102,
                                                       4456, 5, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 6356, 3970, 4102,
                                                       4566, 5, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 6446, 4036, 4102,
                                                       4676, 5, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 6536, 4102, 4786, 5, nmax);

        simdtrf::compute_hrr_geom_100x_dp_out_of_first(buffer, coordinates, 6626, 4168, 4366,
                                                       4896, 5, nmax);

        simdtrf::compute_hrr_geom_100y_dp_out_of_first(buffer, coordinates, 6716, 4234, 4366,
                                                       5006, 5, nmax);

        simdtrf::compute_hrr_geom_100z_dp_out_of_first(buffer, coordinates, 6806, 4300, 4366,
                                                       5116, 5, nmax);

        simdtrf::compute_hrr_dp_out_of_first(buffer, coordinates, 6896, 4366, 5226, 5, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 6986, 4456, 4786,
                                                       5366, 5, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 7136, 4566, 4786,
                                                       5531, 5, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 7286, 4676, 4786,
                                                       5696, 5, nmax);

        simdtrf::compute_hrr_geom_100x_fp_out_of_first(buffer, coordinates, 7436, 4896, 5226,
                                                       5861, 5, nmax);

        simdtrf::compute_hrr_geom_100y_fp_out_of_first(buffer, coordinates, 7586, 5006, 5226,
                                                       6026, 5, nmax);

        simdtrf::compute_hrr_geom_100z_fp_out_of_first(buffer, coordinates, 7736, 5116, 5226,
                                                       6191, 5, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 7886, 6266, 6536,
                                                       6986, 5, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 8066, 6356, 6536,
                                                       7136, 5, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 8246, 6446, 6536,
                                                       7286, 5, nmax);

        simdtrf::compute_hrr_geom_100x_dd_out_of_first(buffer, coordinates, 8426, 6626, 6896,
                                                       7436, 5, nmax);

        simdtrf::compute_hrr_geom_100y_dd_out_of_first(buffer, coordinates, 8606, 6716, 6896,
                                                       7586, 5, nmax);

        simdtrf::compute_hrr_geom_100z_dd_out_of_first(buffer, coordinates, 8786, 6806, 6896,
                                                       7736, 5, nmax);

        simdtrf::transform_d_inner(buffer, 8966, 8426, 6, 5, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 8966, 25, nmax);

        simdtrf::transform_d_inner(buffer, 8966, 8606, 6, 5, nmax);

        simdtrf::transform_d_outer(values + 125 * nvalues + n * npairs, nvalues, buffer, 8966,
                                   25, nmax);

        simdtrf::transform_d_inner(buffer, 8966, 8786, 6, 5, nmax);

        simdtrf::transform_d_outer(values + 250 * nvalues + n * npairs, nvalues, buffer, 8966,
                                   25, nmax);

        simdtrf::transform_d_inner(buffer, 8966, 7886, 6, 5, nmax);

        simdtrf::transform_d_outer(values + 375 * nvalues + n * npairs, nvalues, buffer, 8966,
                                   25, nmax);

        simdtrf::transform_d_inner(buffer, 8966, 8066, 6, 5, nmax);

        simdtrf::transform_d_outer(values + 500 * nvalues + n * npairs, nvalues, buffer, 8966,
                                   25, nmax);

        simdtrf::transform_d_inner(buffer, 8966, 8246, 6, 5, nmax);

        simdtrf::transform_d_outer(values + 625 * nvalues + n * npairs, nvalues, buffer, 8966,
                                   25, nmax);
    }

    for (size_t m = 0; m < 750; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
