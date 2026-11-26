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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecFSG.hpp"

#include <algorithm>
#include <cstddef>
#include <string>

#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "ScreeningFunc.hpp"
#include "SimdDimensions.hpp"
#include "SimdPrimitives.hpp"
#include "SimdBoysFunc.hpp"

#include "SimdGeometryS1.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_fsg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_fsg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 16621, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 16621, 15631, 900, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto beta = b_exps[j];

                const auto fovl = a_norms[i] * b_norms[j];

                const auto fb = a_exps[i] / p;

                const auto fa = -b_exps[j] / p;

                const auto fc = b_exps[j] / p;

                simdfunc::compute_pa(buffer, coordinates, 0, nmax, fa);

                simdfunc::compute_pb(buffer, coordinates, 3, nmax, fb);

                simdfunc::compute_pc(buffer, coordinates, c_coordinates, 6, n, nmax, fc);

                for (size_t k = 0; k < nprim_c; k++)
                {
                    const auto ncols = dimensions[(i * nprim_b + j) * nprim_c + k];

                    if (ncols == 0) continue;

                    const auto gamma = c_exps[k];

                    const auto q = p + gamma;

                    const auto fq = p * gamma / q;

                    const auto fj = 2.0 * fovl * c_norms[k] * pi * pi * std::sqrt(pi)
                                    / (p * gamma * std::sqrt(q));

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 9, 6, 8,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 19, 6, 8,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 3, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 56, 3, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 59, 3, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 62, 3, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 65, 3, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 68, 3, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 71, 3, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 74, 3, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 77, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 80, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 83, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 86, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 89, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 92, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 95, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 98, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 101, 0, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 104, 0, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 107, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 110, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 113, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 116, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 119, 0, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 122, 0, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 125, 0, 6, 10, 11,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 134, 0, 6, 11, 12,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 143, 0, 6, 12, 13,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 152, 0, 6, 13, 14,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 161, 0, 6, 14, 15,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 170, 0, 6, 15, 16,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 179, 0, 6, 16, 17,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 188, 0, 6, 20, 21,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 197, 0, 6, 21, 22,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 206, 0, 6, 22, 23,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 215, 0, 6, 23, 24,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 224, 0, 6, 24, 25,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 233, 0, 6, 25, 26,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 242, 0, 6, 26, 27,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 251, 0, 6, 10, 11,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 257, 0, 6, 11, 12,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 263, 0, 6, 12, 13,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 269, 0, 6, 13, 14,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 275, 0, 6, 14, 15,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 281, 0, 6, 15, 16,
                                                                       92, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 287, 0, 6, 16, 17,
                                                                       95, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 293, 0, 6, 20, 21,
                                                                       101, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 299, 0, 6, 21, 22,
                                                                       104, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 305, 0, 6, 22, 23,
                                                                       107, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 311, 0, 6, 23, 24,
                                                                       110, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 317, 0, 6, 24, 25,
                                                                       113, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 323, 0, 6, 25, 26,
                                                                       116, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 329, 0, 6, 26, 27,
                                                                       119, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 335, 0, 3, 6, 77,
                                                                       80, 125, 134, 251, 257,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 353, 0, 3, 6, 80,
                                                                       83, 134, 143, 257, 263,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 371, 0, 3, 6, 83,
                                                                       86, 143, 152, 263, 269,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 389, 0, 3, 6, 86,
                                                                       89, 152, 161, 269, 275,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 407, 0, 3, 6, 89,
                                                                       92, 161, 170, 275, 281,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 425, 0, 3, 6, 92,
                                                                       95, 170, 179, 281, 287,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 443, 0, 3, 6, 101,
                                                                       104, 188, 197, 293, 299,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 461, 0, 3, 6, 104,
                                                                       107, 197, 206, 299, 305,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 479, 0, 3, 6, 107,
                                                                       110, 206, 215, 305, 311,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 497, 0, 3, 6, 110,
                                                                       113, 215, 224, 311, 317,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 515, 0, 3, 6, 113,
                                                                       116, 224, 233, 317, 323,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 533, 0, 3, 6, 116,
                                                                       119, 233, 242, 323, 329,
                                                                       ncols, gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 551, 0, 6, 77, 80,
                                                                       251, 257, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 561, 0, 6, 80, 83,
                                                                       257, 263, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 571, 0, 6, 83, 86,
                                                                       263, 269, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 581, 0, 6, 86, 89,
                                                                       269, 275, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 591, 0, 6, 89, 92,
                                                                       275, 281, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 601, 0, 6, 92, 95,
                                                                       281, 287, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 611, 0, 6, 101,
                                                                       104, 293, 299, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 621, 0, 6, 104,
                                                                       107, 299, 305, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 631, 0, 6, 107,
                                                                       110, 305, 311, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 641, 0, 6, 110,
                                                                       113, 311, 317, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 651, 0, 6, 113,
                                                                       116, 317, 323, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 661, 0, 6, 116,
                                                                       119, 323, 329, ncols,
                                                                       gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 671, 0, 3, 6, 251,
                                                                       257, 335, 353, 551, 561,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 701, 0, 3, 6, 257,
                                                                       263, 353, 371, 561, 571,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 731, 0, 3, 6, 263,
                                                                       269, 371, 389, 571, 581,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 761, 0, 3, 6, 269,
                                                                       275, 389, 407, 581, 591,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 791, 0, 3, 6, 275,
                                                                       281, 407, 425, 591, 601,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 821, 0, 3, 6, 293,
                                                                       299, 443, 461, 611, 621,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 851, 0, 3, 6, 299,
                                                                       305, 461, 479, 621, 631,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 881, 0, 3, 6, 305,
                                                                       311, 479, 497, 631, 641,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 911, 0, 3, 6, 311,
                                                                       317, 497, 515, 641, 651,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 941, 0, 3, 6, 317,
                                                                       323, 515, 533, 651, 661,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 971, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 974, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 977, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 980, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 983, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 986, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 989, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 992, 6, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 995, 6, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 998, 6, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1001, 6, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1004, 6, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1007, 6, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1010, 6, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1013, 6, 12, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1022, 6, 13, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1031, 6, 14, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1040, 6, 15, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1049, 6, 16, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1058, 6, 17, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1067, 6, 22, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1076, 6, 23, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1085, 6, 24, 65,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1094, 6, 25, 68,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1103, 6, 26, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1112, 6, 27, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1121, 6, 12, 83,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1130, 6, 13, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1139, 6, 14, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1148, 6, 15, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1157, 6, 16, 95,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1166, 6, 17, 98,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1175, 6, 22, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1184, 6, 23, 110,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1193, 6, 24, 113,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1202, 6, 25, 116,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1211, 6, 26, 119,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1220, 6, 27, 122,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1229, 6, 35, 83,
                                                                       143, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1256, 6, 38, 86,
                                                                       152, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1283, 6, 41, 89,
                                                                       161, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1310, 6, 44, 92,
                                                                       170, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1337, 6, 47, 95,
                                                                       179, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1364, 6, 59, 107,
                                                                       206, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1391, 6, 62, 110,
                                                                       215, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1418, 6, 65, 113,
                                                                       224, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1445, 6, 68, 116,
                                                                       233, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1472, 6, 71, 119,
                                                                       242, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1499, 6, 83, 263,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1517, 6, 86, 269,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1535, 6, 89, 275,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1553, 6, 92, 281,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1571, 6, 95, 287,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1589, 6, 107, 305,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1607, 6, 110, 311,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1625, 6, 113, 317,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1643, 6, 116, 323,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1661, 6, 119, 329,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1679, 0, 6, 1229,
                                                                       143, 1256, 263, 371,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1733, 0, 6, 1256,
                                                                       152, 1283, 269, 389,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1787, 0, 6, 1283,
                                                                       161, 1310, 275, 407,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1841, 0, 6, 1310,
                                                                       170, 1337, 281, 425,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1895, 0, 6, 1364,
                                                                       206, 1391, 305, 479,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1949, 0, 6, 1391,
                                                                       215, 1418, 311, 497,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2003, 0, 6, 1418,
                                                                       224, 1445, 317, 515,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2057, 0, 6, 1445,
                                                                       233, 1472, 323, 533,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2111, 6, 263, 571,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2141, 6, 269, 581,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2171, 6, 275, 591,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2201, 6, 281, 601,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2231, 6, 305, 631,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2261, 6, 311, 641,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2291, 6, 317, 651,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2321, 6, 323, 661,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2351, 0, 6, 1679,
                                                                       371, 1733, 571, 731,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2441, 0, 6, 1733,
                                                                       389, 1787, 581, 761,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2531, 0, 6, 1787,
                                                                       407, 1841, 591, 791,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2621, 0, 6, 1895,
                                                                       479, 1949, 631, 881,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2711, 0, 6, 1949,
                                                                       497, 2003, 641, 911,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2801, 0, 6, 2003,
                                                                       515, 2057, 651, 941,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2891, 6, 10, 11,
                                                                       971, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2897, 6, 11, 12,
                                                                       974, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2903, 6, 12, 13,
                                                                       977, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2909, 6, 13, 14,
                                                                       980, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2915, 6, 14, 15,
                                                                       983, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2921, 6, 15, 16,
                                                                       986, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2927, 6, 16, 17,
                                                                       989, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2933, 6, 20, 21,
                                                                       992, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2939, 6, 21, 22,
                                                                       995, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2945, 6, 22, 23,
                                                                       998, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2951, 6, 23, 24,
                                                                       1001, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2957, 6, 24, 25,
                                                                       1004, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2963, 6, 25, 26,
                                                                       1007, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2969, 6, 26, 27,
                                                                       1010, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2975, 3, 6, 2891,
                                                                       971, 2897, 1013, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2993, 3, 6, 2897,
                                                                       974, 2903, 1022, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3011, 3, 6, 2903,
                                                                       977, 2909, 1031, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3029, 3, 6, 2909,
                                                                       980, 2915, 1040, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3047, 3, 6, 2915,
                                                                       983, 2921, 1049, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3065, 3, 6, 2921,
                                                                       986, 2927, 1058, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3083, 3, 6, 2933,
                                                                       992, 2939, 1067, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3101, 3, 6, 2939,
                                                                       995, 2945, 1076, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3119, 3, 6, 2945,
                                                                       998, 2951, 1085, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3137, 3, 6, 2951,
                                                                       1001, 2957, 1094, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3155, 3, 6, 2957,
                                                                       1004, 2963, 1103, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3173, 3, 6, 2963,
                                                                       1007, 2969, 1112, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3191, 0, 6, 2891,
                                                                       971, 2897, 1121, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3209, 0, 6, 2897,
                                                                       974, 2903, 1130, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3227, 0, 6, 2903,
                                                                       977, 2909, 1139, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3245, 0, 6, 2909,
                                                                       980, 2915, 1148, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3263, 0, 6, 2915,
                                                                       983, 2921, 1157, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3281, 0, 6, 2921,
                                                                       986, 2927, 1166, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3299, 0, 6, 2933,
                                                                       992, 2939, 1175, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3317, 0, 6, 2939,
                                                                       995, 2945, 1184, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3335, 0, 6, 2945,
                                                                       998, 2951, 1193, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3353, 0, 6, 2951,
                                                                       1001, 2957, 1202, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3371, 0, 6, 2957,
                                                                       1004, 2963, 1211, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3389, 0, 6, 2963,
                                                                       1007, 2969, 1220, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3407, 0, 3, 6,
                                                                       2975, 1013, 2993, 3191,
                                                                       1121, 3209, 125, 134,
                                                                       1229, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3461, 0, 3, 6,
                                                                       2993, 1022, 3011, 3209,
                                                                       1130, 3227, 134, 143,
                                                                       1256, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3515, 0, 3, 6,
                                                                       3011, 1031, 3029, 3227,
                                                                       1139, 3245, 143, 152,
                                                                       1283, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3569, 0, 3, 6,
                                                                       3029, 1040, 3047, 3245,
                                                                       1148, 3263, 152, 161,
                                                                       1310, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3623, 0, 3, 6,
                                                                       3047, 1049, 3065, 3263,
                                                                       1157, 3281, 161, 170,
                                                                       1337, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3677, 0, 3, 6,
                                                                       3083, 1067, 3101, 3299,
                                                                       1175, 3317, 188, 197,
                                                                       1364, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3731, 0, 3, 6,
                                                                       3101, 1076, 3119, 3317,
                                                                       1184, 3335, 197, 206,
                                                                       1391, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3785, 0, 3, 6,
                                                                       3119, 1085, 3137, 3335,
                                                                       1193, 3353, 206, 215,
                                                                       1418, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3839, 0, 3, 6,
                                                                       3137, 1094, 3155, 3353,
                                                                       1202, 3371, 215, 224,
                                                                       1445, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3893, 0, 3, 6,
                                                                       3155, 1103, 3173, 3371,
                                                                       1211, 3389, 224, 233,
                                                                       1472, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3947, 0, 6, 3191,
                                                                       1121, 3209, 251, 257,
                                                                       1499, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3983, 0, 6, 3209,
                                                                       1130, 3227, 257, 263,
                                                                       1517, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4019, 0, 6, 3227,
                                                                       1139, 3245, 263, 269,
                                                                       1535, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4055, 0, 6, 3245,
                                                                       1148, 3263, 269, 275,
                                                                       1553, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4091, 0, 6, 3263,
                                                                       1157, 3281, 275, 281,
                                                                       1571, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4127, 0, 6, 3299,
                                                                       1175, 3317, 293, 299,
                                                                       1589, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4163, 0, 6, 3317,
                                                                       1184, 3335, 299, 305,
                                                                       1607, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4199, 0, 6, 3335,
                                                                       1193, 3353, 305, 311,
                                                                       1625, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4235, 0, 6, 3353,
                                                                       1202, 3371, 311, 317,
                                                                       1643, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4271, 0, 6, 3371,
                                                                       1211, 3389, 317, 323,
                                                                       1661, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4307, 0, 3, 6,
                                                                       3407, 1229, 3461, 3947,
                                                                       1499, 3983, 335, 353,
                                                                       1679, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4415, 0, 3, 6,
                                                                       3461, 1256, 3515, 3983,
                                                                       1517, 4019, 353, 371,
                                                                       1733, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4523, 0, 3, 6,
                                                                       3515, 1283, 3569, 4019,
                                                                       1535, 4055, 371, 389,
                                                                       1787, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4631, 0, 3, 6,
                                                                       3569, 1310, 3623, 4055,
                                                                       1553, 4091, 389, 407,
                                                                       1841, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4739, 0, 3, 6,
                                                                       3677, 1364, 3731, 4127,
                                                                       1589, 4163, 443, 461,
                                                                       1895, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4847, 0, 3, 6,
                                                                       3731, 1391, 3785, 4163,
                                                                       1607, 4199, 461, 479,
                                                                       1949, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4955, 0, 3, 6,
                                                                       3785, 1418, 3839, 4199,
                                                                       1625, 4235, 479, 497,
                                                                       2003, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5063, 0, 3, 6,
                                                                       3839, 1445, 3893, 4235,
                                                                       1643, 4271, 497, 515,
                                                                       2057, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5171, 0, 6, 3947,
                                                                       1499, 3983, 551, 561,
                                                                       2111, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5231, 0, 6, 3983,
                                                                       1517, 4019, 561, 571,
                                                                       2141, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5291, 0, 6, 4019,
                                                                       1535, 4055, 571, 581,
                                                                       2171, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5351, 0, 6, 4055,
                                                                       1553, 4091, 581, 591,
                                                                       2201, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5411, 0, 6, 4127,
                                                                       1589, 4163, 611, 621,
                                                                       2231, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5471, 0, 6, 4163,
                                                                       1607, 4199, 621, 631,
                                                                       2261, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5531, 0, 6, 4199,
                                                                       1625, 4235, 631, 641,
                                                                       2291, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5591, 0, 6, 4235,
                                                                       1643, 4271, 641, 651,
                                                                       2321, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 5651, 0, 3, 6,
                                                                       3407, 3461, 4307, 1679,
                                                                       4415, 5171, 2111, 5231,
                                                                       671, 701, 2351, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 5831, 0, 3, 6,
                                                                       3461, 3515, 4415, 1733,
                                                                       4523, 5231, 2141, 5291,
                                                                       701, 731, 2441, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6011, 0, 3, 6,
                                                                       3515, 3569, 4523, 1787,
                                                                       4631, 5291, 2171, 5351,
                                                                       731, 761, 2531, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6191, 0, 3, 6,
                                                                       3677, 3731, 4739, 1895,
                                                                       4847, 5411, 2231, 5471,
                                                                       821, 851, 2621, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6371, 0, 3, 6,
                                                                       3731, 3785, 4847, 1949,
                                                                       4955, 5471, 2261, 5531,
                                                                       851, 881, 2711, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6551, 0, 3, 6,
                                                                       3785, 3839, 4955, 2003,
                                                                       5063, 5531, 2291, 5591,
                                                                       881, 911, 2801, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6731, 6, 971, 974,
                                                                       2903, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6741, 6, 974, 977,
                                                                       2909, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6751, 6, 977, 980,
                                                                       2915, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6761, 6, 980, 983,
                                                                       2921, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6771, 6, 983, 986,
                                                                       2927, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6781, 6, 992, 995,
                                                                       2945, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6791, 6, 995, 998,
                                                                       2951, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6801, 6, 998,
                                                                       1001, 2957, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6811, 6, 1001,
                                                                       1004, 2963, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6821, 6, 1004,
                                                                       1007, 2969, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6831, 3, 6, 6731,
                                                                       2903, 6741, 3011, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6861, 3, 6, 6741,
                                                                       2909, 6751, 3029, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6891, 3, 6, 6751,
                                                                       2915, 6761, 3047, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6921, 3, 6, 6761,
                                                                       2921, 6771, 3065, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6951, 3, 6, 6781,
                                                                       2945, 6791, 3119, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6981, 3, 6, 6791,
                                                                       2951, 6801, 3137, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7011, 3, 6, 6801,
                                                                       2957, 6811, 3155, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7041, 3, 6, 6811,
                                                                       2963, 6821, 3173, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7071, 0, 6, 6731,
                                                                       2903, 6741, 1121, 1130,
                                                                       3227, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7101, 0, 6, 6741,
                                                                       2909, 6751, 1130, 1139,
                                                                       3245, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7131, 0, 6, 6751,
                                                                       2915, 6761, 1139, 1148,
                                                                       3263, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7161, 0, 6, 6761,
                                                                       2921, 6771, 1148, 1157,
                                                                       3281, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7191, 0, 6, 6781,
                                                                       2945, 6791, 1175, 1184,
                                                                       3335, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7221, 0, 6, 6791,
                                                                       2951, 6801, 1184, 1193,
                                                                       3353, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7251, 0, 6, 6801,
                                                                       2957, 6811, 1193, 1202,
                                                                       3371, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7281, 0, 6, 6811,
                                                                       2963, 6821, 1202, 1211,
                                                                       3389, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7311, 0, 3, 6,
                                                                       6831, 3011, 6861, 7071,
                                                                       3227, 7101, 1229, 1256,
                                                                       3515, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7401, 0, 3, 6,
                                                                       6861, 3029, 6891, 7101,
                                                                       3245, 7131, 1256, 1283,
                                                                       3569, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7491, 0, 3, 6,
                                                                       6891, 3047, 6921, 7131,
                                                                       3263, 7161, 1283, 1310,
                                                                       3623, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7581, 0, 3, 6,
                                                                       6951, 3119, 6981, 7191,
                                                                       3335, 7221, 1364, 1391,
                                                                       3785, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7671, 0, 3, 6,
                                                                       6981, 3137, 7011, 7221,
                                                                       3353, 7251, 1391, 1418,
                                                                       3839, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 7761, 0, 3, 6,
                                                                       7011, 3155, 7041, 7251,
                                                                       3371, 7281, 1418, 1445,
                                                                       3893, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7851, 0, 6, 7071,
                                                                       3227, 7101, 1499, 1517,
                                                                       4019, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7911, 0, 6, 7101,
                                                                       3245, 7131, 1517, 1535,
                                                                       4055, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 7971, 0, 6, 7131,
                                                                       3263, 7161, 1535, 1553,
                                                                       4091, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8031, 0, 6, 7191,
                                                                       3335, 7221, 1589, 1607,
                                                                       4199, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8091, 0, 6, 7221,
                                                                       3353, 7251, 1607, 1625,
                                                                       4235, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 8151, 0, 6, 7251,
                                                                       3371, 7281, 1625, 1643,
                                                                       4271, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 8211, 0, 3, 6,
                                                                       7311, 3515, 7401, 7851,
                                                                       4019, 7911, 1679, 1733,
                                                                       4523, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 8391, 0, 3, 6,
                                                                       7401, 3569, 7491, 7911,
                                                                       4055, 7971, 1733, 1787,
                                                                       4631, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 8571, 0, 3, 6,
                                                                       7581, 3785, 7671, 8031,
                                                                       4199, 8091, 1895, 1949,
                                                                       4955, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 8751, 0, 3, 6,
                                                                       7671, 3839, 7761, 8091,
                                                                       4235, 8151, 1949, 2003,
                                                                       5063, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 8931, 0, 6, 7851,
                                                                       4019, 7911, 2111, 2141,
                                                                       5291, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9031, 0, 6, 7911,
                                                                       4055, 7971, 2141, 2171,
                                                                       5351, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9131, 0, 6, 8031,
                                                                       4199, 8091, 2231, 2261,
                                                                       5531, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 9231, 0, 6, 8091,
                                                                       4235, 8151, 2261, 2291,
                                                                       5591, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 9331, 0, 3, 6,
                                                                       7311, 7401, 8211, 4523,
                                                                       8391, 8931, 5291, 9031,
                                                                       2351, 2441, 6011, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 9631, 0, 3, 6,
                                                                       7581, 7671, 8571, 4955,
                                                                       8751, 9131, 5531, 9231,
                                                                       2621, 2711, 6551, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9931, 6, 2891,
                                                                       2897, 6731, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9946, 6, 2897,
                                                                       2903, 6741, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9961, 6, 2903,
                                                                       2909, 6751, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9976, 6, 2909,
                                                                       2915, 6761, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9991, 6, 2915,
                                                                       2921, 6771, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10006, 6, 2933,
                                                                       2939, 6781, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10021, 6, 2939,
                                                                       2945, 6791, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10036, 6, 2945,
                                                                       2951, 6801, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10051, 6, 2951,
                                                                       2957, 6811, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10066, 6, 2957,
                                                                       2963, 6821, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10081, 3, 6, 9931,
                                                                       6731, 9946, 2975, 2993,
                                                                       6831, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10126, 3, 6, 9946,
                                                                       6741, 9961, 2993, 3011,
                                                                       6861, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10171, 3, 6, 9961,
                                                                       6751, 9976, 3011, 3029,
                                                                       6891, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10216, 3, 6, 9976,
                                                                       6761, 9991, 3029, 3047,
                                                                       6921, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10261, 3, 6,
                                                                       10006, 6781, 10021, 3083,
                                                                       3101, 6951, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10306, 3, 6,
                                                                       10021, 6791, 10036, 3101,
                                                                       3119, 6981, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10351, 3, 6,
                                                                       10036, 6801, 10051, 3119,
                                                                       3137, 7011, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10396, 3, 6,
                                                                       10051, 6811, 10066, 3137,
                                                                       3155, 7041, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10441, 0, 6, 9931,
                                                                       6731, 9946, 3191, 3209,
                                                                       7071, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10486, 0, 6, 9946,
                                                                       6741, 9961, 3209, 3227,
                                                                       7101, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10531, 0, 6, 9961,
                                                                       6751, 9976, 3227, 3245,
                                                                       7131, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10576, 0, 6, 9976,
                                                                       6761, 9991, 3245, 3263,
                                                                       7161, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10621, 0, 6,
                                                                       10006, 6781, 10021, 3299,
                                                                       3317, 7191, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10666, 0, 6,
                                                                       10021, 6791, 10036, 3317,
                                                                       3335, 7221, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10711, 0, 6,
                                                                       10036, 6801, 10051, 3335,
                                                                       3353, 7251, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 10756, 0, 6,
                                                                       10051, 6811, 10066, 3353,
                                                                       3371, 7281, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 10801, 0, 3, 6,
                                                                       10081, 6831, 10126, 10441,
                                                                       7071, 10486, 3407, 3461,
                                                                       7311, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 10936, 0, 3, 6,
                                                                       10126, 6861, 10171, 10486,
                                                                       7101, 10531, 3461, 3515,
                                                                       7401, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 11071, 0, 3, 6,
                                                                       10171, 6891, 10216, 10531,
                                                                       7131, 10576, 3515, 3569,
                                                                       7491, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 11206, 0, 3, 6,
                                                                       10261, 6951, 10306, 10621,
                                                                       7191, 10666, 3677, 3731,
                                                                       7581, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 11341, 0, 3, 6,
                                                                       10306, 6981, 10351, 10666,
                                                                       7221, 10711, 3731, 3785,
                                                                       7671, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 11476, 0, 3, 6,
                                                                       10351, 7011, 10396, 10711,
                                                                       7251, 10756, 3785, 3839,
                                                                       7761, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11611, 0, 6,
                                                                       10441, 7071, 10486, 3947,
                                                                       3983, 7851, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11701, 0, 6,
                                                                       10486, 7101, 10531, 3983,
                                                                       4019, 7911, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11791, 0, 6,
                                                                       10531, 7131, 10576, 4019,
                                                                       4055, 7971, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11881, 0, 6,
                                                                       10621, 7191, 10666, 4127,
                                                                       4163, 8031, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 11971, 0, 6,
                                                                       10666, 7221, 10711, 4163,
                                                                       4199, 8091, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 12061, 0, 6,
                                                                       10711, 7251, 10756, 4199,
                                                                       4235, 8151, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 12151, 0, 3, 6,
                                                                       10801, 7311, 10936, 11611,
                                                                       7851, 11701, 4307, 4415,
                                                                       8211, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 12421, 0, 3, 6,
                                                                       10936, 7401, 11071, 11701,
                                                                       7911, 11791, 4415, 4523,
                                                                       8391, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 12691, 0, 3, 6,
                                                                       11206, 7581, 11341, 11881,
                                                                       8031, 11971, 4739, 4847,
                                                                       8571, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 12961, 0, 3, 6,
                                                                       11341, 7671, 11476, 11971,
                                                                       8091, 12061, 4847, 4955,
                                                                       8751, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13231, 0, 6,
                                                                       11611, 7851, 11701, 5171,
                                                                       5231, 8931, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13381, 0, 6,
                                                                       11701, 7911, 11791, 5231,
                                                                       5291, 9031, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13531, 0, 6,
                                                                       11881, 8031, 11971, 5411,
                                                                       5471, 9131, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 13681, 0, 6,
                                                                       11971, 8091, 12061, 5471,
                                                                       5531, 9231, ncols, gamma,
                                                                       p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 13831, 0, 3, 6,
                                                                       10801, 10936, 12151, 8211,
                                                                       12421, 13231, 8931, 13381,
                                                                       5651, 5831, 9331, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 14281, 0, 3, 6,
                                                                       11206, 11341, 12691, 8571,
                                                                       12961, 13531, 9131, 13681,
                                                                       6191, 6371, 9631, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 14731, 14281, 10, 15, ncols, beta);

                    simdgeo::geom_s_y(buffer, 14881, 14281, 10, 15, ncols, beta);

                    simdgeo::geom_s_z(buffer, 15031, 14281, 10, 15, ncols, beta);

                    simdgeo::geom_s_x(buffer, 15181, 13831, 10, 15, ncols, beta);

                    simdgeo::geom_s_y(buffer, 15331, 13831, 10, 15, ncols, beta);

                    simdgeo::geom_s_z(buffer, 15481, 13831, 10, 15, ncols, beta);

                    simdfunc::contract_primitives(buffer, 15631, 14731, 900, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 16531, 15631, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 16531, 9, nmax);

        simdtrf::transform_g_inner(buffer, 16531, 15781, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 16531, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 16531, 15931, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 126 * nvalues + n * npairs, nvalues, buffer, 16531,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 16531, 16081, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 189 * nvalues + n * npairs, nvalues, buffer, 16531,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 16531, 16231, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 252 * nvalues + n * npairs, nvalues, buffer, 16531,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 16531, 16381, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 16531,
                                   9, nmax);
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
