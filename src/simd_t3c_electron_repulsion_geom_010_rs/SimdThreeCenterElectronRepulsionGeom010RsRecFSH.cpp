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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecFSH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_fsh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_fsh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 28205, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 462 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 28205, 26835, 1260, dimensions);

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

                    simdfunc::compute_t3c_erf_boys_function(buffer, coordinates, 9, 6, {1, 2, 3,
                                                            4, 5, 6, 7, 8, 9}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 19, 6, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9}, ncols, fj,
                                                        i * nprim_b + j, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 971, 6, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 974, 6, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 977, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 980, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 983, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 986, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 989, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 992, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 995, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 998, 6, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1001, 6, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1004, 6, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1007, 6, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1010, 6, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1013, 6, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1016, 6, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1019, 6, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1022, 6, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1025, 6, 12, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1034, 6, 13, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1043, 6, 14, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1052, 6, 15, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1061, 6, 16, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1070, 6, 17, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1079, 6, 22, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1088, 6, 23, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1097, 6, 24, 65,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1106, 6, 25, 68,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1115, 6, 26, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1124, 6, 27, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1133, 6, 10, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1142, 6, 11, 80,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1151, 6, 12, 83,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1160, 6, 13, 86,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1169, 6, 14, 89,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1178, 6, 15, 92,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1187, 6, 16, 95,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1196, 6, 17, 98,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1205, 6, 20, 101,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1214, 6, 21, 104,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1223, 6, 22, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1232, 6, 23, 110,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1241, 6, 24, 113,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1250, 6, 25, 116,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1259, 6, 26, 119,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1268, 6, 27, 122,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1277, 6, 29, 77,
                                                                       125, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1304, 6, 32, 80,
                                                                       134, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1331, 6, 35, 83,
                                                                       143, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1358, 6, 38, 86,
                                                                       152, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1385, 6, 41, 89,
                                                                       161, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1412, 6, 44, 92,
                                                                       170, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1439, 6, 47, 95,
                                                                       179, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1466, 6, 53, 101,
                                                                       188, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1493, 6, 56, 104,
                                                                       197, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1520, 6, 59, 107,
                                                                       206, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1547, 6, 62, 110,
                                                                       215, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1574, 6, 65, 113,
                                                                       224, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1601, 6, 68, 116,
                                                                       233, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1628, 6, 71, 119,
                                                                       242, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1655, 6, 77, 251,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1673, 6, 80, 257,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1691, 6, 83, 263,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1709, 6, 86, 269,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1727, 6, 89, 275,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1745, 6, 92, 281,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1763, 6, 95, 287,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1781, 6, 101, 293,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1799, 6, 104, 299,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1817, 6, 107, 305,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1835, 6, 110, 311,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1853, 6, 113, 317,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1871, 6, 116, 323,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1889, 6, 119, 329,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1907, 0, 6, 1277,
                                                                       125, 1304, 251, 335,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1961, 0, 6, 1304,
                                                                       134, 1331, 257, 353,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2015, 0, 6, 1331,
                                                                       143, 1358, 263, 371,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2069, 0, 6, 1358,
                                                                       152, 1385, 269, 389,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2123, 0, 6, 1385,
                                                                       161, 1412, 275, 407,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2177, 0, 6, 1412,
                                                                       170, 1439, 281, 425,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2231, 0, 6, 1466,
                                                                       188, 1493, 293, 443,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2285, 0, 6, 1493,
                                                                       197, 1520, 299, 461,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2339, 0, 6, 1520,
                                                                       206, 1547, 305, 479,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2393, 0, 6, 1547,
                                                                       215, 1574, 311, 497,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2447, 0, 6, 1574,
                                                                       224, 1601, 317, 515,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2501, 0, 6, 1601,
                                                                       233, 1628, 323, 533,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2555, 6, 251, 551,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2585, 6, 257, 561,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2615, 6, 263, 571,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2645, 6, 269, 581,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2675, 6, 275, 591,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2705, 6, 281, 601,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2735, 6, 293, 611,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2765, 6, 299, 621,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2795, 6, 305, 631,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2825, 6, 311, 641,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2855, 6, 317, 651,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2885, 6, 323, 661,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2915, 0, 6, 1907,
                                                                       335, 1961, 551, 671,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3005, 0, 6, 1961,
                                                                       353, 2015, 561, 701,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3095, 0, 6, 2015,
                                                                       371, 2069, 571, 731,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3185, 0, 6, 2069,
                                                                       389, 2123, 581, 761,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3275, 0, 6, 2123,
                                                                       407, 2177, 591, 791,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3365, 0, 6, 2231,
                                                                       443, 2285, 611, 821,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3455, 0, 6, 2285,
                                                                       461, 2339, 621, 851,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3545, 0, 6, 2339,
                                                                       479, 2393, 631, 881,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3635, 0, 6, 2393,
                                                                       497, 2447, 641, 911,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3725, 0, 6, 2447,
                                                                       515, 2501, 651, 941,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3815, 6, 10, 11,
                                                                       977, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3821, 6, 11, 12,
                                                                       980, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3827, 6, 12, 13,
                                                                       983, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3833, 6, 13, 14,
                                                                       986, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3839, 6, 14, 15,
                                                                       989, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3845, 6, 15, 16,
                                                                       992, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3851, 6, 16, 17,
                                                                       995, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3857, 6, 20, 21,
                                                                       1004, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3863, 6, 21, 22,
                                                                       1007, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3869, 6, 22, 23,
                                                                       1010, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3875, 6, 23, 24,
                                                                       1013, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3881, 6, 24, 25,
                                                                       1016, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3887, 6, 25, 26,
                                                                       1019, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3893, 6, 26, 27,
                                                                       1022, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3899, 3, 6, 3815,
                                                                       977, 3821, 1025, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3917, 3, 6, 3821,
                                                                       980, 3827, 1034, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3935, 3, 6, 3827,
                                                                       983, 3833, 1043, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3953, 3, 6, 3833,
                                                                       986, 3839, 1052, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3971, 3, 6, 3839,
                                                                       989, 3845, 1061, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3989, 3, 6, 3845,
                                                                       992, 3851, 1070, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4007, 3, 6, 3857,
                                                                       1004, 3863, 1079, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4025, 3, 6, 3863,
                                                                       1007, 3869, 1088, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4043, 3, 6, 3869,
                                                                       1010, 3875, 1097, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4061, 3, 6, 3875,
                                                                       1013, 3881, 1106, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4079, 3, 6, 3881,
                                                                       1016, 3887, 1115, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4097, 3, 6, 3887,
                                                                       1019, 3893, 1124, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4115, 0, 6, 3815,
                                                                       977, 3821, 1151, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4133, 0, 6, 3821,
                                                                       980, 3827, 1160, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4151, 0, 6, 3827,
                                                                       983, 3833, 1169, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4169, 0, 6, 3833,
                                                                       986, 3839, 1178, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4187, 0, 6, 3839,
                                                                       989, 3845, 1187, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4205, 0, 6, 3845,
                                                                       992, 3851, 1196, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4223, 0, 6, 3857,
                                                                       1004, 3863, 1223, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4241, 0, 6, 3863,
                                                                       1007, 3869, 1232, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4259, 0, 6, 3869,
                                                                       1010, 3875, 1241, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4277, 0, 6, 3875,
                                                                       1013, 3881, 1250, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4295, 0, 6, 3881,
                                                                       1016, 3887, 1259, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4313, 0, 6, 3887,
                                                                       1019, 3893, 1268, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4331, 0, 3, 6,
                                                                       3899, 1025, 3917, 4115,
                                                                       1151, 4133, 125, 134,
                                                                       1331, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4385, 0, 3, 6,
                                                                       3917, 1034, 3935, 4133,
                                                                       1160, 4151, 134, 143,
                                                                       1358, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4439, 0, 3, 6,
                                                                       3935, 1043, 3953, 4151,
                                                                       1169, 4169, 143, 152,
                                                                       1385, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4493, 0, 3, 6,
                                                                       3953, 1052, 3971, 4169,
                                                                       1178, 4187, 152, 161,
                                                                       1412, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4547, 0, 3, 6,
                                                                       3971, 1061, 3989, 4187,
                                                                       1187, 4205, 161, 170,
                                                                       1439, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4601, 0, 3, 6,
                                                                       4007, 1079, 4025, 4223,
                                                                       1223, 4241, 188, 197,
                                                                       1520, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4655, 0, 3, 6,
                                                                       4025, 1088, 4043, 4241,
                                                                       1232, 4259, 197, 206,
                                                                       1547, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4709, 0, 3, 6,
                                                                       4043, 1097, 4061, 4259,
                                                                       1241, 4277, 206, 215,
                                                                       1574, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4763, 0, 3, 6,
                                                                       4061, 1106, 4079, 4277,
                                                                       1250, 4295, 215, 224,
                                                                       1601, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4817, 0, 3, 6,
                                                                       4079, 1115, 4097, 4295,
                                                                       1259, 4313, 224, 233,
                                                                       1628, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4871, 0, 6, 4115,
                                                                       1151, 4133, 251, 257,
                                                                       1691, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4907, 0, 6, 4133,
                                                                       1160, 4151, 257, 263,
                                                                       1709, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4943, 0, 6, 4151,
                                                                       1169, 4169, 263, 269,
                                                                       1727, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4979, 0, 6, 4169,
                                                                       1178, 4187, 269, 275,
                                                                       1745, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5015, 0, 6, 4187,
                                                                       1187, 4205, 275, 281,
                                                                       1763, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5051, 0, 6, 4223,
                                                                       1223, 4241, 293, 299,
                                                                       1817, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5087, 0, 6, 4241,
                                                                       1232, 4259, 299, 305,
                                                                       1835, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5123, 0, 6, 4259,
                                                                       1241, 4277, 305, 311,
                                                                       1853, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5159, 0, 6, 4277,
                                                                       1250, 4295, 311, 317,
                                                                       1871, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5195, 0, 6, 4295,
                                                                       1259, 4313, 317, 323,
                                                                       1889, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5231, 0, 3, 6,
                                                                       4331, 1331, 4385, 4871,
                                                                       1691, 4907, 335, 353,
                                                                       2015, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5339, 0, 3, 6,
                                                                       4385, 1358, 4439, 4907,
                                                                       1709, 4943, 353, 371,
                                                                       2069, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5447, 0, 3, 6,
                                                                       4439, 1385, 4493, 4943,
                                                                       1727, 4979, 371, 389,
                                                                       2123, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5555, 0, 3, 6,
                                                                       4493, 1412, 4547, 4979,
                                                                       1745, 5015, 389, 407,
                                                                       2177, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5663, 0, 3, 6,
                                                                       4601, 1520, 4655, 5051,
                                                                       1817, 5087, 443, 461,
                                                                       2339, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5771, 0, 3, 6,
                                                                       4655, 1547, 4709, 5087,
                                                                       1835, 5123, 461, 479,
                                                                       2393, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5879, 0, 3, 6,
                                                                       4709, 1574, 4763, 5123,
                                                                       1853, 5159, 479, 497,
                                                                       2447, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5987, 0, 3, 6,
                                                                       4763, 1601, 4817, 5159,
                                                                       1871, 5195, 497, 515,
                                                                       2501, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6095, 0, 6, 4871,
                                                                       1691, 4907, 551, 561,
                                                                       2615, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6155, 0, 6, 4907,
                                                                       1709, 4943, 561, 571,
                                                                       2645, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6215, 0, 6, 4943,
                                                                       1727, 4979, 571, 581,
                                                                       2675, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6275, 0, 6, 4979,
                                                                       1745, 5015, 581, 591,
                                                                       2705, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6335, 0, 6, 5051,
                                                                       1817, 5087, 611, 621,
                                                                       2795, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6395, 0, 6, 5087,
                                                                       1835, 5123, 621, 631,
                                                                       2825, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6455, 0, 6, 5123,
                                                                       1853, 5159, 631, 641,
                                                                       2855, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6515, 0, 6, 5159,
                                                                       1871, 5195, 641, 651,
                                                                       2885, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6575, 0, 3, 6,
                                                                       4331, 4385, 5231, 2015,
                                                                       5339, 6095, 2615, 6155,
                                                                       671, 701, 3095, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6755, 0, 3, 6,
                                                                       4385, 4439, 5339, 2069,
                                                                       5447, 6155, 2645, 6215,
                                                                       701, 731, 3185, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6935, 0, 3, 6,
                                                                       4439, 4493, 5447, 2123,
                                                                       5555, 6215, 2675, 6275,
                                                                       731, 761, 3275, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 7115, 0, 3, 6,
                                                                       4601, 4655, 5663, 2339,
                                                                       5771, 6335, 2795, 6395,
                                                                       821, 851, 3545, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 7295, 0, 3, 6,
                                                                       4655, 4709, 5771, 2393,
                                                                       5879, 6395, 2825, 6455,
                                                                       851, 881, 3635, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 7475, 0, 3, 6,
                                                                       4709, 4763, 5879, 2447,
                                                                       5987, 6455, 2855, 6515,
                                                                       881, 911, 3725, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7655, 6, 971, 974,
                                                                       3815, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7665, 6, 974, 977,
                                                                       3821, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7675, 6, 977, 980,
                                                                       3827, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7685, 6, 980, 983,
                                                                       3833, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7695, 6, 983, 986,
                                                                       3839, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7705, 6, 986, 989,
                                                                       3845, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7715, 6, 989, 992,
                                                                       3851, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7725, 6, 998,
                                                                       1001, 3857, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7735, 6, 1001,
                                                                       1004, 3863, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7745, 6, 1004,
                                                                       1007, 3869, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7755, 6, 1007,
                                                                       1010, 3875, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7765, 6, 1010,
                                                                       1013, 3881, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7775, 6, 1013,
                                                                       1016, 3887, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 7785, 6, 1016,
                                                                       1019, 3893, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7795, 3, 6, 7655,
                                                                       3815, 7665, 3899, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7825, 3, 6, 7665,
                                                                       3821, 7675, 3917, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7855, 3, 6, 7675,
                                                                       3827, 7685, 3935, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7885, 3, 6, 7685,
                                                                       3833, 7695, 3953, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7915, 3, 6, 7695,
                                                                       3839, 7705, 3971, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7945, 3, 6, 7705,
                                                                       3845, 7715, 3989, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7975, 3, 6, 7725,
                                                                       3857, 7735, 4007, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8005, 3, 6, 7735,
                                                                       3863, 7745, 4025, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8035, 3, 6, 7745,
                                                                       3869, 7755, 4043, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8065, 3, 6, 7755,
                                                                       3875, 7765, 4061, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8095, 3, 6, 7765,
                                                                       3881, 7775, 4079, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8125, 3, 6, 7775,
                                                                       3887, 7785, 4097, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8155, 0, 6, 7655,
                                                                       3815, 7665, 1133, 1142,
                                                                       4115, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8185, 0, 6, 7665,
                                                                       3821, 7675, 1142, 1151,
                                                                       4133, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8215, 0, 6, 7675,
                                                                       3827, 7685, 1151, 1160,
                                                                       4151, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8245, 0, 6, 7685,
                                                                       3833, 7695, 1160, 1169,
                                                                       4169, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8275, 0, 6, 7695,
                                                                       3839, 7705, 1169, 1178,
                                                                       4187, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8305, 0, 6, 7705,
                                                                       3845, 7715, 1178, 1187,
                                                                       4205, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8335, 0, 6, 7725,
                                                                       3857, 7735, 1205, 1214,
                                                                       4223, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8365, 0, 6, 7735,
                                                                       3863, 7745, 1214, 1223,
                                                                       4241, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8395, 0, 6, 7745,
                                                                       3869, 7755, 1223, 1232,
                                                                       4259, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8425, 0, 6, 7755,
                                                                       3875, 7765, 1232, 1241,
                                                                       4277, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8455, 0, 6, 7765,
                                                                       3881, 7775, 1241, 1250,
                                                                       4295, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8485, 0, 6, 7775,
                                                                       3887, 7785, 1250, 1259,
                                                                       4313, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 8515, 0, 3, 6,
                                                                       7795, 3899, 7825, 8155,
                                                                       4115, 8185, 1277, 1304,
                                                                       4331, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 8605, 0, 3, 6,
                                                                       7825, 3917, 7855, 8185,
                                                                       4133, 8215, 1304, 1331,
                                                                       4385, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 8695, 0, 3, 6,
                                                                       7855, 3935, 7885, 8215,
                                                                       4151, 8245, 1331, 1358,
                                                                       4439, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 8785, 0, 3, 6,
                                                                       7885, 3953, 7915, 8245,
                                                                       4169, 8275, 1358, 1385,
                                                                       4493, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 8875, 0, 3, 6,
                                                                       7915, 3971, 7945, 8275,
                                                                       4187, 8305, 1385, 1412,
                                                                       4547, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 8965, 0, 3, 6,
                                                                       7975, 4007, 8005, 8335,
                                                                       4223, 8365, 1466, 1493,
                                                                       4601, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9055, 0, 3, 6,
                                                                       8005, 4025, 8035, 8365,
                                                                       4241, 8395, 1493, 1520,
                                                                       4655, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9145, 0, 3, 6,
                                                                       8035, 4043, 8065, 8395,
                                                                       4259, 8425, 1520, 1547,
                                                                       4709, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9235, 0, 3, 6,
                                                                       8065, 4061, 8095, 8425,
                                                                       4277, 8455, 1547, 1574,
                                                                       4763, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9325, 0, 3, 6,
                                                                       8095, 4079, 8125, 8455,
                                                                       4295, 8485, 1574, 1601,
                                                                       4817, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9415, 0, 6, 8155,
                                                                       4115, 8185, 1655, 1673,
                                                                       4871, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9475, 0, 6, 8185,
                                                                       4133, 8215, 1673, 1691,
                                                                       4907, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9535, 0, 6, 8215,
                                                                       4151, 8245, 1691, 1709,
                                                                       4943, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9595, 0, 6, 8245,
                                                                       4169, 8275, 1709, 1727,
                                                                       4979, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9655, 0, 6, 8275,
                                                                       4187, 8305, 1727, 1745,
                                                                       5015, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9715, 0, 6, 8335,
                                                                       4223, 8365, 1781, 1799,
                                                                       5051, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9775, 0, 6, 8365,
                                                                       4241, 8395, 1799, 1817,
                                                                       5087, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9835, 0, 6, 8395,
                                                                       4259, 8425, 1817, 1835,
                                                                       5123, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9895, 0, 6, 8425,
                                                                       4277, 8455, 1835, 1853,
                                                                       5159, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9955, 0, 6, 8455,
                                                                       4295, 8485, 1853, 1871,
                                                                       5195, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10015, 0, 3, 6,
                                                                       8515, 4331, 8605, 9415,
                                                                       4871, 9475, 1907, 1961,
                                                                       5231, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10195, 0, 3, 6,
                                                                       8605, 4385, 8695, 9475,
                                                                       4907, 9535, 1961, 2015,
                                                                       5339, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10375, 0, 3, 6,
                                                                       8695, 4439, 8785, 9535,
                                                                       4943, 9595, 2015, 2069,
                                                                       5447, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10555, 0, 3, 6,
                                                                       8785, 4493, 8875, 9595,
                                                                       4979, 9655, 2069, 2123,
                                                                       5555, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10735, 0, 3, 6,
                                                                       8965, 4601, 9055, 9715,
                                                                       5051, 9775, 2231, 2285,
                                                                       5663, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10915, 0, 3, 6,
                                                                       9055, 4655, 9145, 9775,
                                                                       5087, 9835, 2285, 2339,
                                                                       5771, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 11095, 0, 3, 6,
                                                                       9145, 4709, 9235, 9835,
                                                                       5123, 9895, 2339, 2393,
                                                                       5879, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 11275, 0, 3, 6,
                                                                       9235, 4763, 9325, 9895,
                                                                       5159, 9955, 2393, 2447,
                                                                       5987, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11455, 0, 6, 9415,
                                                                       4871, 9475, 2555, 2585,
                                                                       6095, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11555, 0, 6, 9475,
                                                                       4907, 9535, 2585, 2615,
                                                                       6155, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11655, 0, 6, 9535,
                                                                       4943, 9595, 2615, 2645,
                                                                       6215, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11755, 0, 6, 9595,
                                                                       4979, 9655, 2645, 2675,
                                                                       6275, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11855, 0, 6, 9715,
                                                                       5051, 9775, 2735, 2765,
                                                                       6335, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11955, 0, 6, 9775,
                                                                       5087, 9835, 2765, 2795,
                                                                       6395, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 12055, 0, 6, 9835,
                                                                       5123, 9895, 2795, 2825,
                                                                       6455, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 12155, 0, 6, 9895,
                                                                       5159, 9955, 2825, 2855,
                                                                       6515, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 12255, 0, 3, 6,
                                                                       8515, 8605, 10015, 5231,
                                                                       10195, 11455, 6095, 11555,
                                                                       2915, 3005, 6575, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 12555, 0, 3, 6,
                                                                       8605, 8695, 10195, 5339,
                                                                       10375, 11555, 6155, 11655,
                                                                       3005, 3095, 6755, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 12855, 0, 3, 6,
                                                                       8695, 8785, 10375, 5447,
                                                                       10555, 11655, 6215, 11755,
                                                                       3095, 3185, 6935, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 13155, 0, 3, 6,
                                                                       8965, 9055, 10735, 5663,
                                                                       10915, 11855, 6335, 11955,
                                                                       3365, 3455, 7115, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 13455, 0, 3, 6,
                                                                       9055, 9145, 10915, 5771,
                                                                       11095, 11955, 6395, 12055,
                                                                       3455, 3545, 7295, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 13755, 0, 3, 6,
                                                                       9145, 9235, 11095, 5879,
                                                                       11275, 12055, 6455, 12155,
                                                                       3545, 3635, 7475, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14055, 6, 3815,
                                                                       3821, 7675, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14070, 6, 3821,
                                                                       3827, 7685, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14085, 6, 3827,
                                                                       3833, 7695, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14100, 6, 3833,
                                                                       3839, 7705, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14115, 6, 3839,
                                                                       3845, 7715, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14130, 6, 3857,
                                                                       3863, 7745, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14145, 6, 3863,
                                                                       3869, 7755, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14160, 6, 3869,
                                                                       3875, 7765, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14175, 6, 3875,
                                                                       3881, 7775, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14190, 6, 3881,
                                                                       3887, 7785, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 14205, 3, 6,
                                                                       14055, 7675, 14070, 3899,
                                                                       3917, 7855, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 14250, 3, 6,
                                                                       14070, 7685, 14085, 3917,
                                                                       3935, 7885, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 14295, 3, 6,
                                                                       14085, 7695, 14100, 3935,
                                                                       3953, 7915, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 14340, 3, 6,
                                                                       14100, 7705, 14115, 3953,
                                                                       3971, 7945, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 14385, 3, 6,
                                                                       14130, 7745, 14145, 4007,
                                                                       4025, 8035, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 14430, 3, 6,
                                                                       14145, 7755, 14160, 4025,
                                                                       4043, 8065, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 14475, 3, 6,
                                                                       14160, 7765, 14175, 4043,
                                                                       4061, 8095, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 14520, 3, 6,
                                                                       14175, 7775, 14190, 4061,
                                                                       4079, 8125, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14565, 0, 6,
                                                                       14055, 7675, 14070, 4115,
                                                                       4133, 8215, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14610, 0, 6,
                                                                       14070, 7685, 14085, 4133,
                                                                       4151, 8245, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14655, 0, 6,
                                                                       14085, 7695, 14100, 4151,
                                                                       4169, 8275, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14700, 0, 6,
                                                                       14100, 7705, 14115, 4169,
                                                                       4187, 8305, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14745, 0, 6,
                                                                       14130, 7745, 14145, 4223,
                                                                       4241, 8395, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14790, 0, 6,
                                                                       14145, 7755, 14160, 4241,
                                                                       4259, 8425, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14835, 0, 6,
                                                                       14160, 7765, 14175, 4259,
                                                                       4277, 8455, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 14880, 0, 6,
                                                                       14175, 7775, 14190, 4277,
                                                                       4295, 8485, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 14925, 0, 3, 6,
                                                                       14205, 7855, 14250, 14565,
                                                                       8215, 14610, 4331, 4385,
                                                                       8695, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 15060, 0, 3, 6,
                                                                       14250, 7885, 14295, 14610,
                                                                       8245, 14655, 4385, 4439,
                                                                       8785, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 15195, 0, 3, 6,
                                                                       14295, 7915, 14340, 14655,
                                                                       8275, 14700, 4439, 4493,
                                                                       8875, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 15330, 0, 3, 6,
                                                                       14385, 8035, 14430, 14745,
                                                                       8395, 14790, 4601, 4655,
                                                                       9145, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 15465, 0, 3, 6,
                                                                       14430, 8065, 14475, 14790,
                                                                       8425, 14835, 4655, 4709,
                                                                       9235, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 15600, 0, 3, 6,
                                                                       14475, 8095, 14520, 14835,
                                                                       8455, 14880, 4709, 4763,
                                                                       9325, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 15735, 0, 6,
                                                                       14565, 8215, 14610, 4871,
                                                                       4907, 9535, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 15825, 0, 6,
                                                                       14610, 8245, 14655, 4907,
                                                                       4943, 9595, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 15915, 0, 6,
                                                                       14655, 8275, 14700, 4943,
                                                                       4979, 9655, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 16005, 0, 6,
                                                                       14745, 8395, 14790, 5051,
                                                                       5087, 9835, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 16095, 0, 6,
                                                                       14790, 8425, 14835, 5087,
                                                                       5123, 9895, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 16185, 0, 6,
                                                                       14835, 8455, 14880, 5123,
                                                                       5159, 9955, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 16275, 0, 3, 6,
                                                                       14925, 8695, 15060, 15735,
                                                                       9535, 15825, 5231, 5339,
                                                                       10375, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 16545, 0, 3, 6,
                                                                       15060, 8785, 15195, 15825,
                                                                       9595, 15915, 5339, 5447,
                                                                       10555, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 16815, 0, 3, 6,
                                                                       15330, 9145, 15465, 16005,
                                                                       9835, 16095, 5663, 5771,
                                                                       11095, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 17085, 0, 3, 6,
                                                                       15465, 9235, 15600, 16095,
                                                                       9895, 16185, 5771, 5879,
                                                                       11275, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 17355, 0, 6,
                                                                       15735, 9535, 15825, 6095,
                                                                       6155, 11655, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 17505, 0, 6,
                                                                       15825, 9595, 15915, 6155,
                                                                       6215, 11755, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 17655, 0, 6,
                                                                       16005, 9835, 16095, 6335,
                                                                       6395, 12055, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 17805, 0, 6,
                                                                       16095, 9895, 16185, 6395,
                                                                       6455, 12155, ncols, gamma,
                                                                       p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 17955, 0, 3, 6,
                                                                       14925, 15060, 16275,
                                                                       10375, 16545, 17355,
                                                                       11655, 17505, 6575, 6755,
                                                                       12855, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 18405, 0, 3, 6,
                                                                       15330, 15465, 16815,
                                                                       11095, 17085, 17655,
                                                                       12055, 17805, 7115, 7295,
                                                                       13755, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 18855, 6, 7655,
                                                                       7665, 14055, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 18876, 6, 7665,
                                                                       7675, 14070, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 18897, 6, 7675,
                                                                       7685, 14085, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 18918, 6, 7685,
                                                                       7695, 14100, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 18939, 6, 7695,
                                                                       7705, 14115, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 18960, 6, 7725,
                                                                       7735, 14130, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 18981, 6, 7735,
                                                                       7745, 14145, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 19002, 6, 7745,
                                                                       7755, 14160, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 19023, 6, 7755,
                                                                       7765, 14175, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 19044, 6, 7765,
                                                                       7775, 14190, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 19065, 3, 6,
                                                                       18855, 14055, 18876, 7795,
                                                                       7825, 14205, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 19128, 3, 6,
                                                                       18876, 14070, 18897, 7825,
                                                                       7855, 14250, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 19191, 3, 6,
                                                                       18897, 14085, 18918, 7855,
                                                                       7885, 14295, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 19254, 3, 6,
                                                                       18918, 14100, 18939, 7885,
                                                                       7915, 14340, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 19317, 3, 6,
                                                                       18960, 14130, 18981, 7975,
                                                                       8005, 14385, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 19380, 3, 6,
                                                                       18981, 14145, 19002, 8005,
                                                                       8035, 14430, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 19443, 3, 6,
                                                                       19002, 14160, 19023, 8035,
                                                                       8065, 14475, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 19506, 3, 6,
                                                                       19023, 14175, 19044, 8065,
                                                                       8095, 14520, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 19569, 0, 6,
                                                                       18855, 14055, 18876, 8155,
                                                                       8185, 14565, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 19632, 0, 6,
                                                                       18876, 14070, 18897, 8185,
                                                                       8215, 14610, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 19695, 0, 6,
                                                                       18897, 14085, 18918, 8215,
                                                                       8245, 14655, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 19758, 0, 6,
                                                                       18918, 14100, 18939, 8245,
                                                                       8275, 14700, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 19821, 0, 6,
                                                                       18960, 14130, 18981, 8335,
                                                                       8365, 14745, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 19884, 0, 6,
                                                                       18981, 14145, 19002, 8365,
                                                                       8395, 14790, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 19947, 0, 6,
                                                                       19002, 14160, 19023, 8395,
                                                                       8425, 14835, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 20010, 0, 6,
                                                                       19023, 14175, 19044, 8425,
                                                                       8455, 14880, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 20073, 0, 3, 6,
                                                                       19065, 14205, 19128,
                                                                       19569, 14565, 19632, 8515,
                                                                       8605, 14925, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 20262, 0, 3, 6,
                                                                       19128, 14250, 19191,
                                                                       19632, 14610, 19695, 8605,
                                                                       8695, 15060, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 20451, 0, 3, 6,
                                                                       19191, 14295, 19254,
                                                                       19695, 14655, 19758, 8695,
                                                                       8785, 15195, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 20640, 0, 3, 6,
                                                                       19317, 14385, 19380,
                                                                       19821, 14745, 19884, 8965,
                                                                       9055, 15330, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 20829, 0, 3, 6,
                                                                       19380, 14430, 19443,
                                                                       19884, 14790, 19947, 9055,
                                                                       9145, 15465, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 21018, 0, 3, 6,
                                                                       19443, 14475, 19506,
                                                                       19947, 14835, 20010, 9145,
                                                                       9235, 15600, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 21207, 0, 6,
                                                                       19569, 14565, 19632, 9415,
                                                                       9475, 15735, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 21333, 0, 6,
                                                                       19632, 14610, 19695, 9475,
                                                                       9535, 15825, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 21459, 0, 6,
                                                                       19695, 14655, 19758, 9535,
                                                                       9595, 15915, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 21585, 0, 6,
                                                                       19821, 14745, 19884, 9715,
                                                                       9775, 16005, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 21711, 0, 6,
                                                                       19884, 14790, 19947, 9775,
                                                                       9835, 16095, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 21837, 0, 6,
                                                                       19947, 14835, 20010, 9835,
                                                                       9895, 16185, ncols, gamma,
                                                                       p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 21963, 0, 3, 6,
                                                                       20073, 14925, 20262,
                                                                       21207, 15735, 21333,
                                                                       10015, 10195, 16275,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 22341, 0, 3, 6,
                                                                       20262, 15060, 20451,
                                                                       21333, 15825, 21459,
                                                                       10195, 10375, 16545,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 22719, 0, 3, 6,
                                                                       20640, 15330, 20829,
                                                                       21585, 16005, 21711,
                                                                       10735, 10915, 16815,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 23097, 0, 3, 6,
                                                                       20829, 15465, 21018,
                                                                       21711, 16095, 21837,
                                                                       10915, 11095, 17085,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 23475, 0, 6,
                                                                       21207, 15735, 21333,
                                                                       11455, 11555, 17355,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 23685, 0, 6,
                                                                       21333, 15825, 21459,
                                                                       11555, 11655, 17505,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 23895, 0, 6,
                                                                       21585, 16005, 21711,
                                                                       11855, 11955, 17655,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 24105, 0, 6,
                                                                       21711, 16095, 21837,
                                                                       11955, 12055, 17805,
                                                                       ncols, gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 24315, 0, 3, 6,
                                                                       20073, 20262, 21963,
                                                                       16275, 22341, 23475,
                                                                       17355, 23685, 12255,
                                                                       12555, 17955, ncols,
                                                                       gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 24945, 0, 3, 6,
                                                                       20640, 20829, 22719,
                                                                       16815, 23097, 23895,
                                                                       17655, 24105, 13155,
                                                                       13455, 18405, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 25575, 24945, 10, 21, ncols, beta);

                    simdgeo::geom_s_y(buffer, 25785, 24945, 10, 21, ncols, beta);

                    simdgeo::geom_s_z(buffer, 25995, 24945, 10, 21, ncols, beta);

                    simdgeo::geom_s_x(buffer, 26205, 24315, 10, 21, ncols, beta);

                    simdgeo::geom_s_y(buffer, 26415, 24315, 10, 21, ncols, beta);

                    simdgeo::geom_s_z(buffer, 26625, 24315, 10, 21, ncols, beta);

                    simdfunc::contract_primitives(buffer, 26835, 25575, 1260, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 28095, 26835, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 28095, 11, nmax);

        simdtrf::transform_h_inner(buffer, 28095, 27045, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 77 * nvalues + n * npairs, nvalues, buffer, 28095,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 28095, 27255, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 154 * nvalues + n * npairs, nvalues, buffer, 28095,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 28095, 27465, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 231 * nvalues + n * npairs, nvalues, buffer, 28095,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 28095, 27675, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 308 * nvalues + n * npairs, nvalues, buffer, 28095,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 28095, 27885, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 385 * nvalues + n * npairs, nvalues, buffer, 28095,
                                   11, nmax);
    }

    for (size_t m = 0; m < 462; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
