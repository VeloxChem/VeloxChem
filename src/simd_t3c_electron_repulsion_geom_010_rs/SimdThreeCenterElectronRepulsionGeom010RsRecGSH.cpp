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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGSH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_gsh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_gsh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 52640, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 52640, 50585, 1890, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 20, 6, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 3, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 61, 3, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 64, 3, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 67, 3, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 70, 3, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 73, 3, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 76, 3, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 79, 3, 6, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 82, 3, 6, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 85, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 88, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 91, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 94, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 97, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 100, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 103, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 106, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 109, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 112, 0, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 115, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 118, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 121, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 124, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 127, 0, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 130, 0, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 133, 0, 6, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 136, 0, 6, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 139, 0, 6, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 148, 0, 6, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 157, 0, 6, 12, 13,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 166, 0, 6, 13, 14,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 175, 0, 6, 14, 15,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 184, 0, 6, 15, 16,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 193, 0, 6, 16, 17,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 202, 0, 6, 17, 18,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 211, 0, 6, 21, 22,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 220, 0, 6, 22, 23,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 229, 0, 6, 23, 24,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 238, 0, 6, 24, 25,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 247, 0, 6, 25, 26,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 256, 0, 6, 26, 27,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 265, 0, 6, 27, 28,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 274, 0, 6, 28, 29,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 283, 0, 6, 10, 11,
                                                                       85, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 289, 0, 6, 11, 12,
                                                                       88, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 295, 0, 6, 12, 13,
                                                                       91, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 301, 0, 6, 13, 14,
                                                                       94, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 307, 0, 6, 14, 15,
                                                                       97, 100, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 313, 0, 6, 15, 16,
                                                                       100, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 319, 0, 6, 16, 17,
                                                                       103, 106, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 325, 0, 6, 17, 18,
                                                                       106, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 331, 0, 6, 21, 22,
                                                                       112, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 337, 0, 6, 22, 23,
                                                                       115, 118, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 343, 0, 6, 23, 24,
                                                                       118, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 349, 0, 6, 24, 25,
                                                                       121, 124, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 355, 0, 6, 25, 26,
                                                                       124, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 361, 0, 6, 26, 27,
                                                                       127, 130, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 367, 0, 6, 27, 28,
                                                                       130, 133, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 373, 0, 6, 28, 29,
                                                                       133, 136, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 379, 0, 3, 6, 85,
                                                                       88, 139, 148, 283, 289,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 397, 0, 3, 6, 88,
                                                                       91, 148, 157, 289, 295,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 415, 0, 3, 6, 91,
                                                                       94, 157, 166, 295, 301,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 433, 0, 3, 6, 94,
                                                                       97, 166, 175, 301, 307,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 451, 0, 3, 6, 97,
                                                                       100, 175, 184, 307, 313,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 469, 0, 3, 6, 100,
                                                                       103, 184, 193, 313, 319,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 487, 0, 3, 6, 103,
                                                                       106, 193, 202, 319, 325,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 505, 0, 3, 6, 112,
                                                                       115, 211, 220, 331, 337,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 523, 0, 3, 6, 115,
                                                                       118, 220, 229, 337, 343,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 541, 0, 3, 6, 118,
                                                                       121, 229, 238, 343, 349,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 559, 0, 3, 6, 121,
                                                                       124, 238, 247, 349, 355,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 577, 0, 3, 6, 124,
                                                                       127, 247, 256, 355, 361,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 595, 0, 3, 6, 127,
                                                                       130, 256, 265, 361, 367,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 613, 0, 3, 6, 130,
                                                                       133, 265, 274, 367, 373,
                                                                       ncols, gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 631, 0, 6, 85, 88,
                                                                       283, 289, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 641, 0, 6, 88, 91,
                                                                       289, 295, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 651, 0, 6, 91, 94,
                                                                       295, 301, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 661, 0, 6, 94, 97,
                                                                       301, 307, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 671, 0, 6, 97,
                                                                       100, 307, 313, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 681, 0, 6, 100,
                                                                       103, 313, 319, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 691, 0, 6, 103,
                                                                       106, 319, 325, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 701, 0, 6, 112,
                                                                       115, 331, 337, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 711, 0, 6, 115,
                                                                       118, 337, 343, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 721, 0, 6, 118,
                                                                       121, 343, 349, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 731, 0, 6, 121,
                                                                       124, 349, 355, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 741, 0, 6, 124,
                                                                       127, 355, 361, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 751, 0, 6, 127,
                                                                       130, 361, 367, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 761, 0, 6, 130,
                                                                       133, 367, 373, ncols,
                                                                       gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 771, 0, 3, 6, 283,
                                                                       289, 379, 397, 631, 641,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 801, 0, 3, 6, 289,
                                                                       295, 397, 415, 641, 651,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 831, 0, 3, 6, 295,
                                                                       301, 415, 433, 651, 661,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 861, 0, 3, 6, 301,
                                                                       307, 433, 451, 661, 671,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 891, 0, 3, 6, 307,
                                                                       313, 451, 469, 671, 681,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 921, 0, 3, 6, 313,
                                                                       319, 469, 487, 681, 691,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 951, 0, 3, 6, 331,
                                                                       337, 505, 523, 701, 711,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 981, 0, 3, 6, 337,
                                                                       343, 523, 541, 711, 721,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1011, 0, 3, 6,
                                                                       343, 349, 541, 559, 721,
                                                                       731, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1041, 0, 3, 6,
                                                                       349, 355, 559, 577, 731,
                                                                       741, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1071, 0, 3, 6,
                                                                       355, 361, 577, 595, 741,
                                                                       751, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1101, 0, 3, 6,
                                                                       361, 367, 595, 613, 751,
                                                                       761, ncols, gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1131, 0, 6, 283,
                                                                       289, 631, 641, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1146, 0, 6, 289,
                                                                       295, 641, 651, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1161, 0, 6, 295,
                                                                       301, 651, 661, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1176, 0, 6, 301,
                                                                       307, 661, 671, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1191, 0, 6, 307,
                                                                       313, 671, 681, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1206, 0, 6, 313,
                                                                       319, 681, 691, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1221, 0, 6, 331,
                                                                       337, 701, 711, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1236, 0, 6, 337,
                                                                       343, 711, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1251, 0, 6, 343,
                                                                       349, 721, 731, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1266, 0, 6, 349,
                                                                       355, 731, 741, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1281, 0, 6, 355,
                                                                       361, 741, 751, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1296, 0, 6, 361,
                                                                       367, 751, 761, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1311, 0, 3, 6,
                                                                       379, 397, 631, 641, 771,
                                                                       801, 1131, 1146, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1356, 0, 3, 6,
                                                                       397, 415, 641, 651, 801,
                                                                       831, 1146, 1161, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1401, 0, 3, 6,
                                                                       415, 433, 651, 661, 831,
                                                                       861, 1161, 1176, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1446, 0, 3, 6,
                                                                       433, 451, 661, 671, 861,
                                                                       891, 1176, 1191, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1491, 0, 3, 6,
                                                                       451, 469, 671, 681, 891,
                                                                       921, 1191, 1206, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1536, 0, 3, 6,
                                                                       505, 523, 701, 711, 951,
                                                                       981, 1221, 1236, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1581, 0, 3, 6,
                                                                       523, 541, 711, 721, 981,
                                                                       1011, 1236, 1251, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1626, 0, 3, 6,
                                                                       541, 559, 721, 731, 1011,
                                                                       1041, 1251, 1266, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1671, 0, 3, 6,
                                                                       559, 577, 731, 741, 1041,
                                                                       1071, 1266, 1281, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 6,
                                                                       577, 595, 741, 751, 1071,
                                                                       1101, 1281, 1296, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1761, 6, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1764, 6, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1767, 6, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1770, 6, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1773, 6, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1776, 6, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1779, 6, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1782, 6, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1785, 6, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1788, 6, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1791, 6, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1794, 6, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1797, 6, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1800, 6, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1803, 6, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1806, 6, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1809, 6, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1812, 6, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1815, 6, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1818, 6, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1821, 6, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1830, 6, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1839, 6, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1848, 6, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1857, 6, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1866, 6, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1875, 6, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1884, 6, 23, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1893, 6, 24, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1902, 6, 25, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1911, 6, 26, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1920, 6, 27, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1929, 6, 28, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1938, 6, 29, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1947, 6, 10, 85,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1956, 6, 11, 88,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1965, 6, 12, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1974, 6, 13, 94,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1983, 6, 14, 97,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1992, 6, 15, 100,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2001, 6, 16, 103,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2010, 6, 17, 106,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2019, 6, 18, 109,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2028, 6, 21, 112,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2037, 6, 22, 115,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2046, 6, 23, 118,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2055, 6, 24, 121,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2064, 6, 25, 124,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2073, 6, 26, 127,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2082, 6, 27, 130,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2091, 6, 28, 133,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2100, 6, 29, 136,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2109, 6, 31, 85,
                                                                       139, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2136, 6, 34, 88,
                                                                       148, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2163, 6, 37, 91,
                                                                       157, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2190, 6, 40, 94,
                                                                       166, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2217, 6, 43, 97,
                                                                       175, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2244, 6, 46, 100,
                                                                       184, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2271, 6, 49, 103,
                                                                       193, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2298, 6, 52, 106,
                                                                       202, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2325, 6, 58, 112,
                                                                       211, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2352, 6, 61, 115,
                                                                       220, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2379, 6, 64, 118,
                                                                       229, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2406, 6, 67, 121,
                                                                       238, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2433, 6, 70, 124,
                                                                       247, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2460, 6, 73, 127,
                                                                       256, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2487, 6, 76, 130,
                                                                       265, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2514, 6, 79, 133,
                                                                       274, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2541, 6, 85, 283,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2559, 6, 88, 289,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2577, 6, 91, 295,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2595, 6, 94, 301,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2613, 6, 97, 307,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2631, 6, 100, 313,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2649, 6, 103, 319,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2667, 6, 106, 325,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2685, 6, 112, 331,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2703, 6, 115, 337,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2721, 6, 118, 343,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2739, 6, 121, 349,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2757, 6, 124, 355,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2775, 6, 127, 361,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2793, 6, 130, 367,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2811, 6, 133, 373,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2829, 0, 6, 2109,
                                                                       139, 2136, 283, 379,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2883, 0, 6, 2136,
                                                                       148, 2163, 289, 397,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2937, 0, 6, 2163,
                                                                       157, 2190, 295, 415,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2991, 0, 6, 2190,
                                                                       166, 2217, 301, 433,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3045, 0, 6, 2217,
                                                                       175, 2244, 307, 451,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3099, 0, 6, 2244,
                                                                       184, 2271, 313, 469,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3153, 0, 6, 2271,
                                                                       193, 2298, 319, 487,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3207, 0, 6, 2325,
                                                                       211, 2352, 331, 505,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3261, 0, 6, 2352,
                                                                       220, 2379, 337, 523,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3315, 0, 6, 2379,
                                                                       229, 2406, 343, 541,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3369, 0, 6, 2406,
                                                                       238, 2433, 349, 559,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3423, 0, 6, 2433,
                                                                       247, 2460, 355, 577,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3477, 0, 6, 2460,
                                                                       256, 2487, 361, 595,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3531, 0, 6, 2487,
                                                                       265, 2514, 367, 613,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3585, 6, 283, 631,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3615, 6, 289, 641,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3645, 6, 295, 651,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3675, 6, 301, 661,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3705, 6, 307, 671,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3735, 6, 313, 681,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3765, 6, 319, 691,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3795, 6, 331, 701,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3825, 6, 337, 711,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3855, 6, 343, 721,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3885, 6, 349, 731,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3915, 6, 355, 741,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3945, 6, 361, 751,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3975, 6, 367, 761,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4005, 0, 6, 2829,
                                                                       379, 2883, 631, 771,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4095, 0, 6, 2883,
                                                                       397, 2937, 641, 801,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4185, 0, 6, 2937,
                                                                       415, 2991, 651, 831,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4275, 0, 6, 2991,
                                                                       433, 3045, 661, 861,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4365, 0, 6, 3045,
                                                                       451, 3099, 671, 891,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4455, 0, 6, 3099,
                                                                       469, 3153, 681, 921,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4545, 0, 6, 3207,
                                                                       505, 3261, 701, 951,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4635, 0, 6, 3261,
                                                                       523, 3315, 711, 981,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4725, 0, 6, 3315,
                                                                       541, 3369, 721, 1011,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4815, 0, 6, 3369,
                                                                       559, 3423, 731, 1041,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4905, 0, 6, 3423,
                                                                       577, 3477, 741, 1071,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4995, 0, 6, 3477,
                                                                       595, 3531, 751, 1101,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5085, 6, 631,
                                                                       1131, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5130, 6, 641,
                                                                       1146, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5175, 6, 651,
                                                                       1161, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5220, 6, 661,
                                                                       1176, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5265, 6, 671,
                                                                       1191, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5310, 6, 681,
                                                                       1206, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5355, 6, 701,
                                                                       1221, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5400, 6, 711,
                                                                       1236, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5445, 6, 721,
                                                                       1251, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5490, 6, 731,
                                                                       1266, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5535, 6, 741,
                                                                       1281, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5580, 6, 751,
                                                                       1296, ncols, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 5625, 0, 6, 4005,
                                                                       771, 4095, 1131, 1311,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 5760, 0, 6, 4095,
                                                                       801, 4185, 1146, 1356,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 5895, 0, 6, 4185,
                                                                       831, 4275, 1161, 1401,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 6030, 0, 6, 4275,
                                                                       861, 4365, 1176, 1446,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 6165, 0, 6, 4365,
                                                                       891, 4455, 1191, 1491,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 6300, 0, 6, 4545,
                                                                       951, 4635, 1221, 1536,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 6435, 0, 6, 4635,
                                                                       981, 4725, 1236, 1581,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 6570, 0, 6, 4725,
                                                                       1011, 4815, 1251, 1626,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 6705, 0, 6, 4815,
                                                                       1041, 4905, 1266, 1671,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 6840, 0, 6, 4905,
                                                                       1071, 4995, 1281, 1716,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6975, 6, 10, 11,
                                                                       1767, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6981, 6, 11, 12,
                                                                       1770, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6987, 6, 12, 13,
                                                                       1773, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6993, 6, 13, 14,
                                                                       1776, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6999, 6, 14, 15,
                                                                       1779, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7005, 6, 15, 16,
                                                                       1782, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7011, 6, 16, 17,
                                                                       1785, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7017, 6, 17, 18,
                                                                       1788, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7023, 6, 21, 22,
                                                                       1797, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7029, 6, 22, 23,
                                                                       1800, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7035, 6, 23, 24,
                                                                       1803, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7041, 6, 24, 25,
                                                                       1806, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7047, 6, 25, 26,
                                                                       1809, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7053, 6, 26, 27,
                                                                       1812, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7059, 6, 27, 28,
                                                                       1815, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7065, 6, 28, 29,
                                                                       1818, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7071, 3, 6, 6975,
                                                                       1767, 6981, 1821, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7089, 3, 6, 6981,
                                                                       1770, 6987, 1830, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7107, 3, 6, 6987,
                                                                       1773, 6993, 1839, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7125, 3, 6, 6993,
                                                                       1776, 6999, 1848, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7143, 3, 6, 6999,
                                                                       1779, 7005, 1857, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7161, 3, 6, 7005,
                                                                       1782, 7011, 1866, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7179, 3, 6, 7011,
                                                                       1785, 7017, 1875, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7197, 3, 6, 7023,
                                                                       1797, 7029, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7215, 3, 6, 7029,
                                                                       1800, 7035, 1893, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7233, 3, 6, 7035,
                                                                       1803, 7041, 1902, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7251, 3, 6, 7041,
                                                                       1806, 7047, 1911, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7269, 3, 6, 7047,
                                                                       1809, 7053, 1920, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7287, 3, 6, 7053,
                                                                       1812, 7059, 1929, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7305, 3, 6, 7059,
                                                                       1815, 7065, 1938, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7323, 0, 6, 6975,
                                                                       1767, 6981, 1965, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7341, 0, 6, 6981,
                                                                       1770, 6987, 1974, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7359, 0, 6, 6987,
                                                                       1773, 6993, 1983, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7377, 0, 6, 6993,
                                                                       1776, 6999, 1992, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7395, 0, 6, 6999,
                                                                       1779, 7005, 2001, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7413, 0, 6, 7005,
                                                                       1782, 7011, 2010, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7431, 0, 6, 7011,
                                                                       1785, 7017, 2019, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7449, 0, 6, 7023,
                                                                       1797, 7029, 2046, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7467, 0, 6, 7029,
                                                                       1800, 7035, 2055, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7485, 0, 6, 7035,
                                                                       1803, 7041, 2064, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7503, 0, 6, 7041,
                                                                       1806, 7047, 2073, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7521, 0, 6, 7047,
                                                                       1809, 7053, 2082, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7539, 0, 6, 7053,
                                                                       1812, 7059, 2091, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7557, 0, 6, 7059,
                                                                       1815, 7065, 2100, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7575, 0, 3, 6,
                                                                       7071, 1821, 7089, 7323,
                                                                       1965, 7341, 139, 148,
                                                                       2163, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7629, 0, 3, 6,
                                                                       7089, 1830, 7107, 7341,
                                                                       1974, 7359, 148, 157,
                                                                       2190, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7683, 0, 3, 6,
                                                                       7107, 1839, 7125, 7359,
                                                                       1983, 7377, 157, 166,
                                                                       2217, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7737, 0, 3, 6,
                                                                       7125, 1848, 7143, 7377,
                                                                       1992, 7395, 166, 175,
                                                                       2244, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7791, 0, 3, 6,
                                                                       7143, 1857, 7161, 7395,
                                                                       2001, 7413, 175, 184,
                                                                       2271, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7845, 0, 3, 6,
                                                                       7161, 1866, 7179, 7413,
                                                                       2010, 7431, 184, 193,
                                                                       2298, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7899, 0, 3, 6,
                                                                       7197, 1884, 7215, 7449,
                                                                       2046, 7467, 211, 220,
                                                                       2379, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7953, 0, 3, 6,
                                                                       7215, 1893, 7233, 7467,
                                                                       2055, 7485, 220, 229,
                                                                       2406, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8007, 0, 3, 6,
                                                                       7233, 1902, 7251, 7485,
                                                                       2064, 7503, 229, 238,
                                                                       2433, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8061, 0, 3, 6,
                                                                       7251, 1911, 7269, 7503,
                                                                       2073, 7521, 238, 247,
                                                                       2460, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8115, 0, 3, 6,
                                                                       7269, 1920, 7287, 7521,
                                                                       2082, 7539, 247, 256,
                                                                       2487, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8169, 0, 3, 6,
                                                                       7287, 1929, 7305, 7539,
                                                                       2091, 7557, 256, 265,
                                                                       2514, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8223, 0, 6, 7323,
                                                                       1965, 7341, 283, 289,
                                                                       2577, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8259, 0, 6, 7341,
                                                                       1974, 7359, 289, 295,
                                                                       2595, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8295, 0, 6, 7359,
                                                                       1983, 7377, 295, 301,
                                                                       2613, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8331, 0, 6, 7377,
                                                                       1992, 7395, 301, 307,
                                                                       2631, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8367, 0, 6, 7395,
                                                                       2001, 7413, 307, 313,
                                                                       2649, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8403, 0, 6, 7413,
                                                                       2010, 7431, 313, 319,
                                                                       2667, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8439, 0, 6, 7449,
                                                                       2046, 7467, 331, 337,
                                                                       2721, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8475, 0, 6, 7467,
                                                                       2055, 7485, 337, 343,
                                                                       2739, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8511, 0, 6, 7485,
                                                                       2064, 7503, 343, 349,
                                                                       2757, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8547, 0, 6, 7503,
                                                                       2073, 7521, 349, 355,
                                                                       2775, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8583, 0, 6, 7521,
                                                                       2082, 7539, 355, 361,
                                                                       2793, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 8619, 0, 6, 7539,
                                                                       2091, 7557, 361, 367,
                                                                       2811, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 8655, 0, 3, 6,
                                                                       7575, 2163, 7629, 8223,
                                                                       2577, 8259, 379, 397,
                                                                       2937, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 8763, 0, 3, 6,
                                                                       7629, 2190, 7683, 8259,
                                                                       2595, 8295, 397, 415,
                                                                       2991, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 8871, 0, 3, 6,
                                                                       7683, 2217, 7737, 8295,
                                                                       2613, 8331, 415, 433,
                                                                       3045, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 8979, 0, 3, 6,
                                                                       7737, 2244, 7791, 8331,
                                                                       2631, 8367, 433, 451,
                                                                       3099, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 9087, 0, 3, 6,
                                                                       7791, 2271, 7845, 8367,
                                                                       2649, 8403, 451, 469,
                                                                       3153, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 9195, 0, 3, 6,
                                                                       7899, 2379, 7953, 8439,
                                                                       2721, 8475, 505, 523,
                                                                       3315, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 9303, 0, 3, 6,
                                                                       7953, 2406, 8007, 8475,
                                                                       2739, 8511, 523, 541,
                                                                       3369, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 9411, 0, 3, 6,
                                                                       8007, 2433, 8061, 8511,
                                                                       2757, 8547, 541, 559,
                                                                       3423, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 9519, 0, 3, 6,
                                                                       8061, 2460, 8115, 8547,
                                                                       2775, 8583, 559, 577,
                                                                       3477, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 9627, 0, 3, 6,
                                                                       8115, 2487, 8169, 8583,
                                                                       2793, 8619, 577, 595,
                                                                       3531, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9735, 0, 6, 8223,
                                                                       2577, 8259, 631, 641,
                                                                       3645, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9795, 0, 6, 8259,
                                                                       2595, 8295, 641, 651,
                                                                       3675, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9855, 0, 6, 8295,
                                                                       2613, 8331, 651, 661,
                                                                       3705, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9915, 0, 6, 8331,
                                                                       2631, 8367, 661, 671,
                                                                       3735, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 9975, 0, 6, 8367,
                                                                       2649, 8403, 671, 681,
                                                                       3765, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10035, 0, 6, 8439,
                                                                       2721, 8475, 701, 711,
                                                                       3855, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10095, 0, 6, 8475,
                                                                       2739, 8511, 711, 721,
                                                                       3885, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10155, 0, 6, 8511,
                                                                       2757, 8547, 721, 731,
                                                                       3915, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10215, 0, 6, 8547,
                                                                       2775, 8583, 731, 741,
                                                                       3945, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 10275, 0, 6, 8583,
                                                                       2793, 8619, 741, 751,
                                                                       3975, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 10335, 0, 3, 6,
                                                                       7575, 7629, 8655, 2937,
                                                                       8763, 9735, 3645, 9795,
                                                                       771, 801, 4185, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 10515, 0, 3, 6,
                                                                       7629, 7683, 8763, 2991,
                                                                       8871, 9795, 3675, 9855,
                                                                       801, 831, 4275, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 10695, 0, 3, 6,
                                                                       7683, 7737, 8871, 3045,
                                                                       8979, 9855, 3705, 9915,
                                                                       831, 861, 4365, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 10875, 0, 3, 6,
                                                                       7737, 7791, 8979, 3099,
                                                                       9087, 9915, 3735, 9975,
                                                                       861, 891, 4455, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 11055, 0, 3, 6,
                                                                       7899, 7953, 9195, 3315,
                                                                       9303, 10035, 3855, 10095,
                                                                       951, 981, 4725, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 11235, 0, 3, 6,
                                                                       7953, 8007, 9303, 3369,
                                                                       9411, 10095, 3885, 10155,
                                                                       981, 1011, 4815, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 11415, 0, 3, 6,
                                                                       8007, 8061, 9411, 3423,
                                                                       9519, 10155, 3915, 10215,
                                                                       1011, 1041, 4905, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 11595, 0, 3, 6,
                                                                       8061, 8115, 9519, 3477,
                                                                       9627, 10215, 3945, 10275,
                                                                       1041, 1071, 4995, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11775, 0, 6, 9735,
                                                                       3645, 9795, 1131, 1146,
                                                                       5175, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11865, 0, 6, 9795,
                                                                       3675, 9855, 1146, 1161,
                                                                       5220, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 11955, 0, 6, 9855,
                                                                       3705, 9915, 1161, 1176,
                                                                       5265, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12045, 0, 6, 9915,
                                                                       3735, 9975, 1176, 1191,
                                                                       5310, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12135, 0, 6,
                                                                       10035, 3855, 10095, 1221,
                                                                       1236, 5445, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12225, 0, 6,
                                                                       10095, 3885, 10155, 1236,
                                                                       1251, 5490, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12315, 0, 6,
                                                                       10155, 3915, 10215, 1251,
                                                                       1266, 5535, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 12405, 0, 6,
                                                                       10215, 3945, 10275, 1266,
                                                                       1281, 5580, ncols, gamma,
                                                                       p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 12495, 0, 3, 6,
                                                                       8655, 8763, 10335, 4185,
                                                                       10515, 11775, 5175, 11865,
                                                                       1311, 1356, 5895, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 12765, 0, 3, 6,
                                                                       8763, 8871, 10515, 4275,
                                                                       10695, 11865, 5220, 11955,
                                                                       1356, 1401, 6030, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 13035, 0, 3, 6,
                                                                       8871, 8979, 10695, 4365,
                                                                       10875, 11955, 5265, 12045,
                                                                       1401, 1446, 6165, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 13305, 0, 3, 6,
                                                                       9195, 9303, 11055, 4725,
                                                                       11235, 12135, 5445, 12225,
                                                                       1536, 1581, 6570, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 13575, 0, 3, 6,
                                                                       9303, 9411, 11235, 4815,
                                                                       11415, 12225, 5490, 12315,
                                                                       1581, 1626, 6705, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 13845, 0, 3, 6,
                                                                       9411, 9519, 11415, 4905,
                                                                       11595, 12315, 5535, 12405,
                                                                       1626, 1671, 6840, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14115, 6, 1761,
                                                                       1764, 6975, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14125, 6, 1764,
                                                                       1767, 6981, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14135, 6, 1767,
                                                                       1770, 6987, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14145, 6, 1770,
                                                                       1773, 6993, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14155, 6, 1773,
                                                                       1776, 6999, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14165, 6, 1776,
                                                                       1779, 7005, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14175, 6, 1779,
                                                                       1782, 7011, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14185, 6, 1782,
                                                                       1785, 7017, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14195, 6, 1791,
                                                                       1794, 7023, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14205, 6, 1794,
                                                                       1797, 7029, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14215, 6, 1797,
                                                                       1800, 7035, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14225, 6, 1800,
                                                                       1803, 7041, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14235, 6, 1803,
                                                                       1806, 7047, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14245, 6, 1806,
                                                                       1809, 7053, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14255, 6, 1809,
                                                                       1812, 7059, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14265, 6, 1812,
                                                                       1815, 7065, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14275, 3, 6,
                                                                       14115, 6975, 14125, 7071,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14305, 3, 6,
                                                                       14125, 6981, 14135, 7089,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14335, 3, 6,
                                                                       14135, 6987, 14145, 7107,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14365, 3, 6,
                                                                       14145, 6993, 14155, 7125,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14395, 3, 6,
                                                                       14155, 6999, 14165, 7143,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14425, 3, 6,
                                                                       14165, 7005, 14175, 7161,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14455, 3, 6,
                                                                       14175, 7011, 14185, 7179,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14485, 3, 6,
                                                                       14195, 7023, 14205, 7197,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14515, 3, 6,
                                                                       14205, 7029, 14215, 7215,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14545, 3, 6,
                                                                       14215, 7035, 14225, 7233,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14575, 3, 6,
                                                                       14225, 7041, 14235, 7251,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14605, 3, 6,
                                                                       14235, 7047, 14245, 7269,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14635, 3, 6,
                                                                       14245, 7053, 14255, 7287,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14665, 3, 6,
                                                                       14255, 7059, 14265, 7305,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14695, 0, 6,
                                                                       14115, 6975, 14125, 1947,
                                                                       1956, 7323, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14725, 0, 6,
                                                                       14125, 6981, 14135, 1956,
                                                                       1965, 7341, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14755, 0, 6,
                                                                       14135, 6987, 14145, 1965,
                                                                       1974, 7359, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14785, 0, 6,
                                                                       14145, 6993, 14155, 1974,
                                                                       1983, 7377, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14815, 0, 6,
                                                                       14155, 6999, 14165, 1983,
                                                                       1992, 7395, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14845, 0, 6,
                                                                       14165, 7005, 14175, 1992,
                                                                       2001, 7413, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14875, 0, 6,
                                                                       14175, 7011, 14185, 2001,
                                                                       2010, 7431, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14905, 0, 6,
                                                                       14195, 7023, 14205, 2028,
                                                                       2037, 7449, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14935, 0, 6,
                                                                       14205, 7029, 14215, 2037,
                                                                       2046, 7467, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14965, 0, 6,
                                                                       14215, 7035, 14225, 2046,
                                                                       2055, 7485, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14995, 0, 6,
                                                                       14225, 7041, 14235, 2055,
                                                                       2064, 7503, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15025, 0, 6,
                                                                       14235, 7047, 14245, 2064,
                                                                       2073, 7521, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15055, 0, 6,
                                                                       14245, 7053, 14255, 2073,
                                                                       2082, 7539, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 15085, 0, 6,
                                                                       14255, 7059, 14265, 2082,
                                                                       2091, 7557, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15115, 0, 3, 6,
                                                                       14275, 7071, 14305, 14695,
                                                                       7323, 14725, 2109, 2136,
                                                                       7575, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15205, 0, 3, 6,
                                                                       14305, 7089, 14335, 14725,
                                                                       7341, 14755, 2136, 2163,
                                                                       7629, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15295, 0, 3, 6,
                                                                       14335, 7107, 14365, 14755,
                                                                       7359, 14785, 2163, 2190,
                                                                       7683, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15385, 0, 3, 6,
                                                                       14365, 7125, 14395, 14785,
                                                                       7377, 14815, 2190, 2217,
                                                                       7737, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15475, 0, 3, 6,
                                                                       14395, 7143, 14425, 14815,
                                                                       7395, 14845, 2217, 2244,
                                                                       7791, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15565, 0, 3, 6,
                                                                       14425, 7161, 14455, 14845,
                                                                       7413, 14875, 2244, 2271,
                                                                       7845, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15655, 0, 3, 6,
                                                                       14485, 7197, 14515, 14905,
                                                                       7449, 14935, 2325, 2352,
                                                                       7899, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15745, 0, 3, 6,
                                                                       14515, 7215, 14545, 14935,
                                                                       7467, 14965, 2352, 2379,
                                                                       7953, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15835, 0, 3, 6,
                                                                       14545, 7233, 14575, 14965,
                                                                       7485, 14995, 2379, 2406,
                                                                       8007, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15925, 0, 3, 6,
                                                                       14575, 7251, 14605, 14995,
                                                                       7503, 15025, 2406, 2433,
                                                                       8061, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 16015, 0, 3, 6,
                                                                       14605, 7269, 14635, 15025,
                                                                       7521, 15055, 2433, 2460,
                                                                       8115, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 16105, 0, 3, 6,
                                                                       14635, 7287, 14665, 15055,
                                                                       7539, 15085, 2460, 2487,
                                                                       8169, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16195, 0, 6,
                                                                       14695, 7323, 14725, 2541,
                                                                       2559, 8223, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16255, 0, 6,
                                                                       14725, 7341, 14755, 2559,
                                                                       2577, 8259, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16315, 0, 6,
                                                                       14755, 7359, 14785, 2577,
                                                                       2595, 8295, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16375, 0, 6,
                                                                       14785, 7377, 14815, 2595,
                                                                       2613, 8331, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16435, 0, 6,
                                                                       14815, 7395, 14845, 2613,
                                                                       2631, 8367, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16495, 0, 6,
                                                                       14845, 7413, 14875, 2631,
                                                                       2649, 8403, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16555, 0, 6,
                                                                       14905, 7449, 14935, 2685,
                                                                       2703, 8439, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16615, 0, 6,
                                                                       14935, 7467, 14965, 2703,
                                                                       2721, 8475, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16675, 0, 6,
                                                                       14965, 7485, 14995, 2721,
                                                                       2739, 8511, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16735, 0, 6,
                                                                       14995, 7503, 15025, 2739,
                                                                       2757, 8547, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16795, 0, 6,
                                                                       15025, 7521, 15055, 2757,
                                                                       2775, 8583, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 16855, 0, 6,
                                                                       15055, 7539, 15085, 2775,
                                                                       2793, 8619, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 16915, 0, 3, 6,
                                                                       15115, 7575, 15205, 16195,
                                                                       8223, 16255, 2829, 2883,
                                                                       8655, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 17095, 0, 3, 6,
                                                                       15205, 7629, 15295, 16255,
                                                                       8259, 16315, 2883, 2937,
                                                                       8763, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 17275, 0, 3, 6,
                                                                       15295, 7683, 15385, 16315,
                                                                       8295, 16375, 2937, 2991,
                                                                       8871, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 17455, 0, 3, 6,
                                                                       15385, 7737, 15475, 16375,
                                                                       8331, 16435, 2991, 3045,
                                                                       8979, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 17635, 0, 3, 6,
                                                                       15475, 7791, 15565, 16435,
                                                                       8367, 16495, 3045, 3099,
                                                                       9087, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 17815, 0, 3, 6,
                                                                       15655, 7899, 15745, 16555,
                                                                       8439, 16615, 3207, 3261,
                                                                       9195, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 17995, 0, 3, 6,
                                                                       15745, 7953, 15835, 16615,
                                                                       8475, 16675, 3261, 3315,
                                                                       9303, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 18175, 0, 3, 6,
                                                                       15835, 8007, 15925, 16675,
                                                                       8511, 16735, 3315, 3369,
                                                                       9411, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 18355, 0, 3, 6,
                                                                       15925, 8061, 16015, 16735,
                                                                       8547, 16795, 3369, 3423,
                                                                       9519, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 18535, 0, 3, 6,
                                                                       16015, 8115, 16105, 16795,
                                                                       8583, 16855, 3423, 3477,
                                                                       9627, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18715, 0, 6,
                                                                       16195, 8223, 16255, 3585,
                                                                       3615, 9735, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18815, 0, 6,
                                                                       16255, 8259, 16315, 3615,
                                                                       3645, 9795, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 18915, 0, 6,
                                                                       16315, 8295, 16375, 3645,
                                                                       3675, 9855, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19015, 0, 6,
                                                                       16375, 8331, 16435, 3675,
                                                                       3705, 9915, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19115, 0, 6,
                                                                       16435, 8367, 16495, 3705,
                                                                       3735, 9975, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19215, 0, 6,
                                                                       16555, 8439, 16615, 3795,
                                                                       3825, 10035, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19315, 0, 6,
                                                                       16615, 8475, 16675, 3825,
                                                                       3855, 10095, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19415, 0, 6,
                                                                       16675, 8511, 16735, 3855,
                                                                       3885, 10155, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19515, 0, 6,
                                                                       16735, 8547, 16795, 3885,
                                                                       3915, 10215, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 19615, 0, 6,
                                                                       16795, 8583, 16855, 3915,
                                                                       3945, 10275, ncols, gamma,
                                                                       p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 19715, 0, 3, 6,
                                                                       15115, 15205, 16915, 8655,
                                                                       17095, 18715, 9735, 18815,
                                                                       4005, 4095, 10335, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 20015, 0, 3, 6,
                                                                       15205, 15295, 17095, 8763,
                                                                       17275, 18815, 9795, 18915,
                                                                       4095, 4185, 10515, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 20315, 0, 3, 6,
                                                                       15295, 15385, 17275, 8871,
                                                                       17455, 18915, 9855, 19015,
                                                                       4185, 4275, 10695, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 20615, 0, 3, 6,
                                                                       15385, 15475, 17455, 8979,
                                                                       17635, 19015, 9915, 19115,
                                                                       4275, 4365, 10875, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 20915, 0, 3, 6,
                                                                       15655, 15745, 17815, 9195,
                                                                       17995, 19215, 10035,
                                                                       19315, 4545, 4635, 11055,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 21215, 0, 3, 6,
                                                                       15745, 15835, 17995, 9303,
                                                                       18175, 19315, 10095,
                                                                       19415, 4635, 4725, 11235,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 21515, 0, 3, 6,
                                                                       15835, 15925, 18175, 9411,
                                                                       18355, 19415, 10155,
                                                                       19515, 4725, 4815, 11415,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 21815, 0, 3, 6,
                                                                       15925, 16015, 18355, 9519,
                                                                       18535, 19515, 10215,
                                                                       19615, 4815, 4905, 11595,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22115, 0, 6,
                                                                       18715, 9735, 18815, 5085,
                                                                       5130, 11775, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22265, 0, 6,
                                                                       18815, 9795, 18915, 5130,
                                                                       5175, 11865, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22415, 0, 6,
                                                                       18915, 9855, 19015, 5175,
                                                                       5220, 11955, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22565, 0, 6,
                                                                       19015, 9915, 19115, 5220,
                                                                       5265, 12045, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22715, 0, 6,
                                                                       19215, 10035, 19315, 5355,
                                                                       5400, 12135, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 22865, 0, 6,
                                                                       19315, 10095, 19415, 5400,
                                                                       5445, 12225, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23015, 0, 6,
                                                                       19415, 10155, 19515, 5445,
                                                                       5490, 12315, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 23165, 0, 6,
                                                                       19515, 10215, 19615, 5490,
                                                                       5535, 12405, ncols, gamma,
                                                                       p, q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 23315, 0, 3, 6,
                                                                       16915, 17095, 19715,
                                                                       10335, 20015, 22115,
                                                                       11775, 22265, 5625, 5760,
                                                                       12495, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 23765, 0, 3, 6,
                                                                       17095, 17275, 20015,
                                                                       10515, 20315, 22265,
                                                                       11865, 22415, 5760, 5895,
                                                                       12765, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 24215, 0, 3, 6,
                                                                       17275, 17455, 20315,
                                                                       10695, 20615, 22415,
                                                                       11955, 22565, 5895, 6030,
                                                                       13035, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 24665, 0, 3, 6,
                                                                       17815, 17995, 20915,
                                                                       11055, 21215, 22715,
                                                                       12135, 22865, 6300, 6435,
                                                                       13305, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 25115, 0, 3, 6,
                                                                       17995, 18175, 21215,
                                                                       11235, 21515, 22865,
                                                                       12225, 23015, 6435, 6570,
                                                                       13575, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 25565, 0, 3, 6,
                                                                       18175, 18355, 21515,
                                                                       11415, 21815, 23015,
                                                                       12315, 23165, 6570, 6705,
                                                                       13845, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26015, 6, 6975,
                                                                       6981, 14135, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26030, 6, 6981,
                                                                       6987, 14145, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26045, 6, 6987,
                                                                       6993, 14155, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26060, 6, 6993,
                                                                       6999, 14165, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26075, 6, 6999,
                                                                       7005, 14175, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26090, 6, 7005,
                                                                       7011, 14185, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26105, 6, 7023,
                                                                       7029, 14215, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26120, 6, 7029,
                                                                       7035, 14225, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26135, 6, 7035,
                                                                       7041, 14235, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26150, 6, 7041,
                                                                       7047, 14245, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26165, 6, 7047,
                                                                       7053, 14255, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26180, 6, 7053,
                                                                       7059, 14265, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26195, 3, 6,
                                                                       26015, 14135, 26030, 7071,
                                                                       7089, 14335, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26240, 3, 6,
                                                                       26030, 14145, 26045, 7089,
                                                                       7107, 14365, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26285, 3, 6,
                                                                       26045, 14155, 26060, 7107,
                                                                       7125, 14395, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26330, 3, 6,
                                                                       26060, 14165, 26075, 7125,
                                                                       7143, 14425, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26375, 3, 6,
                                                                       26075, 14175, 26090, 7143,
                                                                       7161, 14455, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26420, 3, 6,
                                                                       26105, 14215, 26120, 7197,
                                                                       7215, 14545, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26465, 3, 6,
                                                                       26120, 14225, 26135, 7215,
                                                                       7233, 14575, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26510, 3, 6,
                                                                       26135, 14235, 26150, 7233,
                                                                       7251, 14605, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26555, 3, 6,
                                                                       26150, 14245, 26165, 7251,
                                                                       7269, 14635, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26600, 3, 6,
                                                                       26165, 14255, 26180, 7269,
                                                                       7287, 14665, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26645, 0, 6,
                                                                       26015, 14135, 26030, 7323,
                                                                       7341, 14755, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26690, 0, 6,
                                                                       26030, 14145, 26045, 7341,
                                                                       7359, 14785, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26735, 0, 6,
                                                                       26045, 14155, 26060, 7359,
                                                                       7377, 14815, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26780, 0, 6,
                                                                       26060, 14165, 26075, 7377,
                                                                       7395, 14845, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26825, 0, 6,
                                                                       26075, 14175, 26090, 7395,
                                                                       7413, 14875, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26870, 0, 6,
                                                                       26105, 14215, 26120, 7449,
                                                                       7467, 14965, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26915, 0, 6,
                                                                       26120, 14225, 26135, 7467,
                                                                       7485, 14995, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 26960, 0, 6,
                                                                       26135, 14235, 26150, 7485,
                                                                       7503, 15025, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27005, 0, 6,
                                                                       26150, 14245, 26165, 7503,
                                                                       7521, 15055, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 27050, 0, 6,
                                                                       26165, 14255, 26180, 7521,
                                                                       7539, 15085, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 27095, 0, 3, 6,
                                                                       26195, 14335, 26240,
                                                                       26645, 14755, 26690, 7575,
                                                                       7629, 15295, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 27230, 0, 3, 6,
                                                                       26240, 14365, 26285,
                                                                       26690, 14785, 26735, 7629,
                                                                       7683, 15385, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 27365, 0, 3, 6,
                                                                       26285, 14395, 26330,
                                                                       26735, 14815, 26780, 7683,
                                                                       7737, 15475, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 27500, 0, 3, 6,
                                                                       26330, 14425, 26375,
                                                                       26780, 14845, 26825, 7737,
                                                                       7791, 15565, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 27635, 0, 3, 6,
                                                                       26420, 14545, 26465,
                                                                       26870, 14965, 26915, 7899,
                                                                       7953, 15835, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 27770, 0, 3, 6,
                                                                       26465, 14575, 26510,
                                                                       26915, 14995, 26960, 7953,
                                                                       8007, 15925, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 27905, 0, 3, 6,
                                                                       26510, 14605, 26555,
                                                                       26960, 15025, 27005, 8007,
                                                                       8061, 16015, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 28040, 0, 3, 6,
                                                                       26555, 14635, 26600,
                                                                       27005, 15055, 27050, 8061,
                                                                       8115, 16105, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28175, 0, 6,
                                                                       26645, 14755, 26690, 8223,
                                                                       8259, 16315, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28265, 0, 6,
                                                                       26690, 14785, 26735, 8259,
                                                                       8295, 16375, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28355, 0, 6,
                                                                       26735, 14815, 26780, 8295,
                                                                       8331, 16435, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28445, 0, 6,
                                                                       26780, 14845, 26825, 8331,
                                                                       8367, 16495, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28535, 0, 6,
                                                                       26870, 14965, 26915, 8439,
                                                                       8475, 16675, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28625, 0, 6,
                                                                       26915, 14995, 26960, 8475,
                                                                       8511, 16735, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28715, 0, 6,
                                                                       26960, 15025, 27005, 8511,
                                                                       8547, 16795, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 28805, 0, 6,
                                                                       27005, 15055, 27050, 8547,
                                                                       8583, 16855, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 28895, 0, 3, 6,
                                                                       27095, 15295, 27230,
                                                                       28175, 16315, 28265, 8655,
                                                                       8763, 17275, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 29165, 0, 3, 6,
                                                                       27230, 15385, 27365,
                                                                       28265, 16375, 28355, 8763,
                                                                       8871, 17455, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 29435, 0, 3, 6,
                                                                       27365, 15475, 27500,
                                                                       28355, 16435, 28445, 8871,
                                                                       8979, 17635, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 29705, 0, 3, 6,
                                                                       27635, 15835, 27770,
                                                                       28535, 16675, 28625, 9195,
                                                                       9303, 18175, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 29975, 0, 3, 6,
                                                                       27770, 15925, 27905,
                                                                       28625, 16735, 28715, 9303,
                                                                       9411, 18355, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 30245, 0, 3, 6,
                                                                       27905, 16015, 28040,
                                                                       28715, 16795, 28805, 9411,
                                                                       9519, 18535, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 30515, 0, 6,
                                                                       28175, 16315, 28265, 9735,
                                                                       9795, 18915, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 30665, 0, 6,
                                                                       28265, 16375, 28355, 9795,
                                                                       9855, 19015, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 30815, 0, 6,
                                                                       28355, 16435, 28445, 9855,
                                                                       9915, 19115, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 30965, 0, 6,
                                                                       28535, 16675, 28625,
                                                                       10035, 10095, 19415,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 31115, 0, 6,
                                                                       28625, 16735, 28715,
                                                                       10095, 10155, 19515,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 31265, 0, 6,
                                                                       28715, 16795, 28805,
                                                                       10155, 10215, 19615,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 31415, 0, 3, 6,
                                                                       27095, 27230, 28895,
                                                                       17275, 29165, 30515,
                                                                       18915, 30665, 10335,
                                                                       10515, 20315, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 31865, 0, 3, 6,
                                                                       27230, 27365, 29165,
                                                                       17455, 29435, 30665,
                                                                       19015, 30815, 10515,
                                                                       10695, 20615, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 32315, 0, 3, 6,
                                                                       27635, 27770, 29705,
                                                                       18175, 29975, 30965,
                                                                       19415, 31115, 11055,
                                                                       11235, 21515, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 32765, 0, 3, 6,
                                                                       27770, 27905, 29975,
                                                                       18355, 30245, 31115,
                                                                       19515, 31265, 11235,
                                                                       11415, 21815, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 33215, 0, 6,
                                                                       30515, 18915, 30665,
                                                                       11775, 11865, 22415,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 33440, 0, 6,
                                                                       30665, 19015, 30815,
                                                                       11865, 11955, 22565,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 33665, 0, 6,
                                                                       30965, 19415, 31115,
                                                                       12135, 12225, 23015,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 33890, 0, 6,
                                                                       31115, 19515, 31265,
                                                                       12225, 12315, 23165,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 34115, 0, 3, 6,
                                                                       28895, 29165, 31415,
                                                                       20315, 31865, 33215,
                                                                       22415, 33440, 12495,
                                                                       12765, 24215, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 34790, 0, 3, 6,
                                                                       29705, 29975, 32315,
                                                                       21515, 32765, 33665,
                                                                       23015, 33890, 13305,
                                                                       13575, 25565, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35465, 6, 14115,
                                                                       14125, 26015, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35486, 6, 14125,
                                                                       14135, 26030, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35507, 6, 14135,
                                                                       14145, 26045, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35528, 6, 14145,
                                                                       14155, 26060, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35549, 6, 14155,
                                                                       14165, 26075, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35570, 6, 14165,
                                                                       14175, 26090, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35591, 6, 14195,
                                                                       14205, 26105, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35612, 6, 14205,
                                                                       14215, 26120, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35633, 6, 14215,
                                                                       14225, 26135, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35654, 6, 14225,
                                                                       14235, 26150, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35675, 6, 14235,
                                                                       14245, 26165, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35696, 6, 14245,
                                                                       14255, 26180, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35717, 3, 6,
                                                                       35465, 26015, 35486,
                                                                       14275, 14305, 26195,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35780, 3, 6,
                                                                       35486, 26030, 35507,
                                                                       14305, 14335, 26240,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35843, 3, 6,
                                                                       35507, 26045, 35528,
                                                                       14335, 14365, 26285,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35906, 3, 6,
                                                                       35528, 26060, 35549,
                                                                       14365, 14395, 26330,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35969, 3, 6,
                                                                       35549, 26075, 35570,
                                                                       14395, 14425, 26375,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36032, 3, 6,
                                                                       35591, 26105, 35612,
                                                                       14485, 14515, 26420,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36095, 3, 6,
                                                                       35612, 26120, 35633,
                                                                       14515, 14545, 26465,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36158, 3, 6,
                                                                       35633, 26135, 35654,
                                                                       14545, 14575, 26510,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36221, 3, 6,
                                                                       35654, 26150, 35675,
                                                                       14575, 14605, 26555,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36284, 3, 6,
                                                                       35675, 26165, 35696,
                                                                       14605, 14635, 26600,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 36347, 0, 6,
                                                                       35465, 26015, 35486,
                                                                       14695, 14725, 26645,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 36410, 0, 6,
                                                                       35486, 26030, 35507,
                                                                       14725, 14755, 26690,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 36473, 0, 6,
                                                                       35507, 26045, 35528,
                                                                       14755, 14785, 26735,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 36536, 0, 6,
                                                                       35528, 26060, 35549,
                                                                       14785, 14815, 26780,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 36599, 0, 6,
                                                                       35549, 26075, 35570,
                                                                       14815, 14845, 26825,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 36662, 0, 6,
                                                                       35591, 26105, 35612,
                                                                       14905, 14935, 26870,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 36725, 0, 6,
                                                                       35612, 26120, 35633,
                                                                       14935, 14965, 26915,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 36788, 0, 6,
                                                                       35633, 26135, 35654,
                                                                       14965, 14995, 26960,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 36851, 0, 6,
                                                                       35654, 26150, 35675,
                                                                       14995, 15025, 27005,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 36914, 0, 6,
                                                                       35675, 26165, 35696,
                                                                       15025, 15055, 27050,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 36977, 0, 3, 6,
                                                                       35717, 26195, 35780,
                                                                       36347, 26645, 36410,
                                                                       15115, 15205, 27095,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 37166, 0, 3, 6,
                                                                       35780, 26240, 35843,
                                                                       36410, 26690, 36473,
                                                                       15205, 15295, 27230,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 37355, 0, 3, 6,
                                                                       35843, 26285, 35906,
                                                                       36473, 26735, 36536,
                                                                       15295, 15385, 27365,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 37544, 0, 3, 6,
                                                                       35906, 26330, 35969,
                                                                       36536, 26780, 36599,
                                                                       15385, 15475, 27500,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 37733, 0, 3, 6,
                                                                       36032, 26420, 36095,
                                                                       36662, 26870, 36725,
                                                                       15655, 15745, 27635,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 37922, 0, 3, 6,
                                                                       36095, 26465, 36158,
                                                                       36725, 26915, 36788,
                                                                       15745, 15835, 27770,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 38111, 0, 3, 6,
                                                                       36158, 26510, 36221,
                                                                       36788, 26960, 36851,
                                                                       15835, 15925, 27905,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 38300, 0, 3, 6,
                                                                       36221, 26555, 36284,
                                                                       36851, 27005, 36914,
                                                                       15925, 16015, 28040,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38489, 0, 6,
                                                                       36347, 26645, 36410,
                                                                       16195, 16255, 28175,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38615, 0, 6,
                                                                       36410, 26690, 36473,
                                                                       16255, 16315, 28265,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38741, 0, 6,
                                                                       36473, 26735, 36536,
                                                                       16315, 16375, 28355,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38867, 0, 6,
                                                                       36536, 26780, 36599,
                                                                       16375, 16435, 28445,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 38993, 0, 6,
                                                                       36662, 26870, 36725,
                                                                       16555, 16615, 28535,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 39119, 0, 6,
                                                                       36725, 26915, 36788,
                                                                       16615, 16675, 28625,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 39245, 0, 6,
                                                                       36788, 26960, 36851,
                                                                       16675, 16735, 28715,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 39371, 0, 6,
                                                                       36851, 27005, 36914,
                                                                       16735, 16795, 28805,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 39497, 0, 3, 6,
                                                                       36977, 27095, 37166,
                                                                       38489, 28175, 38615,
                                                                       16915, 17095, 28895,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 39875, 0, 3, 6,
                                                                       37166, 27230, 37355,
                                                                       38615, 28265, 38741,
                                                                       17095, 17275, 29165,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 40253, 0, 3, 6,
                                                                       37355, 27365, 37544,
                                                                       38741, 28355, 38867,
                                                                       17275, 17455, 29435,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 40631, 0, 3, 6,
                                                                       37733, 27635, 37922,
                                                                       38993, 28535, 39119,
                                                                       17815, 17995, 29705,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 41009, 0, 3, 6,
                                                                       37922, 27770, 38111,
                                                                       39119, 28625, 39245,
                                                                       17995, 18175, 29975,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 41387, 0, 3, 6,
                                                                       38111, 27905, 38300,
                                                                       39245, 28715, 39371,
                                                                       18175, 18355, 30245,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 41765, 0, 6,
                                                                       38489, 28175, 38615,
                                                                       18715, 18815, 30515,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 41975, 0, 6,
                                                                       38615, 28265, 38741,
                                                                       18815, 18915, 30665,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 42185, 0, 6,
                                                                       38741, 28355, 38867,
                                                                       18915, 19015, 30815,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 42395, 0, 6,
                                                                       38993, 28535, 39119,
                                                                       19215, 19315, 30965,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 42605, 0, 6,
                                                                       39119, 28625, 39245,
                                                                       19315, 19415, 31115,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 42815, 0, 6,
                                                                       39245, 28715, 39371,
                                                                       19415, 19515, 31265,
                                                                       ncols, gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 43025, 0, 3, 6,
                                                                       36977, 37166, 39497,
                                                                       28895, 39875, 41765,
                                                                       30515, 41975, 19715,
                                                                       20015, 31415, ncols,
                                                                       gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 43655, 0, 3, 6,
                                                                       37166, 37355, 39875,
                                                                       29165, 40253, 41975,
                                                                       30665, 42185, 20015,
                                                                       20315, 31865, ncols,
                                                                       gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 44285, 0, 3, 6,
                                                                       37733, 37922, 40631,
                                                                       29705, 41009, 42395,
                                                                       30965, 42605, 20915,
                                                                       21215, 32315, ncols,
                                                                       gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 44915, 0, 3, 6,
                                                                       37922, 38111, 41009,
                                                                       29975, 41387, 42605,
                                                                       31115, 42815, 21215,
                                                                       21515, 32765, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 45545, 0, 6,
                                                                       41765, 30515, 41975,
                                                                       22115, 22265, 33215,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 45860, 0, 6,
                                                                       41975, 30665, 42185,
                                                                       22265, 22415, 33440,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 46175, 0, 6,
                                                                       42395, 30965, 42605,
                                                                       22715, 22865, 33665,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 46490, 0, 6,
                                                                       42605, 31115, 42815,
                                                                       22865, 23015, 33890,
                                                                       ncols, gamma, p, q);

                    compute_prim_gph_three_center_electron_repulsion_0(buffer, 46805, 0, 3, 6,
                                                                       39497, 39875, 43025,
                                                                       31415, 43655, 45545,
                                                                       33215, 45860, 23315,
                                                                       23765, 34115, ncols,
                                                                       gamma, p, q);

                    compute_prim_gph_three_center_electron_repulsion_0(buffer, 47750, 0, 3, 6,
                                                                       40631, 41009, 44285,
                                                                       32315, 44915, 46175,
                                                                       33665, 46490, 24665,
                                                                       25115, 34790, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 48695, 47750, 15, 21, ncols, beta);

                    simdgeo::geom_s_y(buffer, 49010, 47750, 15, 21, ncols, beta);

                    simdgeo::geom_s_z(buffer, 49325, 47750, 15, 21, ncols, beta);

                    simdgeo::geom_s_x(buffer, 49640, 46805, 15, 21, ncols, beta);

                    simdgeo::geom_s_y(buffer, 49955, 46805, 15, 21, ncols, beta);

                    simdgeo::geom_s_z(buffer, 50270, 46805, 15, 21, ncols, beta);

                    simdfunc::contract_primitives(buffer, 50585, 48695, 1890, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 52475, 50585, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 52475, 11, nmax);

        simdtrf::transform_h_inner(buffer, 52475, 50900, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 99 * nvalues + n * npairs, nvalues, buffer, 52475,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 52475, 51215, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 198 * nvalues + n * npairs, nvalues, buffer, 52475,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 52475, 51530, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 297 * nvalues + n * npairs, nvalues, buffer, 52475,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 52475, 51845, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 396 * nvalues + n * npairs, nvalues, buffer, 52475,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 52475, 52160, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 495 * nvalues + n * npairs, nvalues, buffer, 52475,
                                   11, nmax);
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
