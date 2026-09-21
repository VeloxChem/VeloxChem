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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGSG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
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
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_gsg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_gsg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 31056, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 486 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 31056, 29571, 1350, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 9, 6, 9,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 20, 6, 9,
                                                             ncols, fj, i * nprim_b + j, fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1761, 6, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1764, 6, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1767, 6, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1770, 6, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1773, 6, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1776, 6, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1779, 6, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1782, 6, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1785, 6, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1788, 6, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1791, 6, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1794, 6, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1797, 6, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1800, 6, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1803, 6, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1806, 6, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1809, 6, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1818, 6, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1827, 6, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1836, 6, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1845, 6, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1854, 6, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1863, 6, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1872, 6, 23, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1881, 6, 24, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1890, 6, 25, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1899, 6, 26, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1908, 6, 27, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1917, 6, 28, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1926, 6, 29, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1935, 6, 12, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1944, 6, 13, 94,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1953, 6, 14, 97,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1962, 6, 15, 100,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1971, 6, 16, 103,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1980, 6, 17, 106,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1989, 6, 18, 109,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1998, 6, 23, 118,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2007, 6, 24, 121,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2016, 6, 25, 124,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2025, 6, 26, 127,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2034, 6, 27, 130,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2043, 6, 28, 133,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2052, 6, 29, 136,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2061, 6, 37, 91,
                                                                       157, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2088, 6, 40, 94,
                                                                       166, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2115, 6, 43, 97,
                                                                       175, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2142, 6, 46, 100,
                                                                       184, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2169, 6, 49, 103,
                                                                       193, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2196, 6, 52, 106,
                                                                       202, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2223, 6, 64, 118,
                                                                       229, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2250, 6, 67, 121,
                                                                       238, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2277, 6, 70, 124,
                                                                       247, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2304, 6, 73, 127,
                                                                       256, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2331, 6, 76, 130,
                                                                       265, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2358, 6, 79, 133,
                                                                       274, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2385, 6, 91, 295,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2403, 6, 94, 301,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2421, 6, 97, 307,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2439, 6, 100, 313,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2457, 6, 103, 319,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2475, 6, 106, 325,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2493, 6, 118, 343,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2511, 6, 121, 349,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2529, 6, 124, 355,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2547, 6, 127, 361,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2565, 6, 130, 367,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2583, 6, 133, 373,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2601, 0, 6, 2061,
                                                                       157, 2088, 295, 415,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2655, 0, 6, 2088,
                                                                       166, 2115, 301, 433,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2709, 0, 6, 2115,
                                                                       175, 2142, 307, 451,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2763, 0, 6, 2142,
                                                                       184, 2169, 313, 469,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2817, 0, 6, 2169,
                                                                       193, 2196, 319, 487,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2871, 0, 6, 2223,
                                                                       229, 2250, 343, 541,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2925, 0, 6, 2250,
                                                                       238, 2277, 349, 559,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2979, 0, 6, 2277,
                                                                       247, 2304, 355, 577,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3033, 0, 6, 2304,
                                                                       256, 2331, 361, 595,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3087, 0, 6, 2331,
                                                                       265, 2358, 367, 613,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3141, 6, 295, 651,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3171, 6, 301, 661,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3201, 6, 307, 671,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3231, 6, 313, 681,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3261, 6, 319, 691,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3291, 6, 343, 721,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3321, 6, 349, 731,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3351, 6, 355, 741,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3381, 6, 361, 751,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3411, 6, 367, 761,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3441, 0, 6, 2601,
                                                                       415, 2655, 651, 831,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3531, 0, 6, 2655,
                                                                       433, 2709, 661, 861,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3621, 0, 6, 2709,
                                                                       451, 2763, 671, 891,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3711, 0, 6, 2763,
                                                                       469, 2817, 681, 921,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3801, 0, 6, 2871,
                                                                       541, 2925, 721, 1011,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3891, 0, 6, 2925,
                                                                       559, 2979, 731, 1041,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3981, 0, 6, 2979,
                                                                       577, 3033, 741, 1071,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4071, 0, 6, 3033,
                                                                       595, 3087, 751, 1101,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4161, 6, 651,
                                                                       1161, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4206, 6, 661,
                                                                       1176, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4251, 6, 671,
                                                                       1191, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4296, 6, 681,
                                                                       1206, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4341, 6, 721,
                                                                       1251, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4386, 6, 731,
                                                                       1266, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4431, 6, 741,
                                                                       1281, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 4476, 6, 751,
                                                                       1296, ncols, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 4521, 0, 6, 3441,
                                                                       831, 3531, 1161, 1401,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 4656, 0, 6, 3531,
                                                                       861, 3621, 1176, 1446,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 4791, 0, 6, 3621,
                                                                       891, 3711, 1191, 1491,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 4926, 0, 6, 3801,
                                                                       1011, 3891, 1251, 1626,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 5061, 0, 6, 3891,
                                                                       1041, 3981, 1266, 1671,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 5196, 0, 6, 3981,
                                                                       1071, 4071, 1281, 1716,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5331, 6, 10, 11,
                                                                       1761, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5337, 6, 11, 12,
                                                                       1764, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5343, 6, 12, 13,
                                                                       1767, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5349, 6, 13, 14,
                                                                       1770, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5355, 6, 14, 15,
                                                                       1773, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5361, 6, 15, 16,
                                                                       1776, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5367, 6, 16, 17,
                                                                       1779, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5373, 6, 17, 18,
                                                                       1782, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5379, 6, 21, 22,
                                                                       1785, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5385, 6, 22, 23,
                                                                       1788, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5391, 6, 23, 24,
                                                                       1791, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5397, 6, 24, 25,
                                                                       1794, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5403, 6, 25, 26,
                                                                       1797, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5409, 6, 26, 27,
                                                                       1800, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5415, 6, 27, 28,
                                                                       1803, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5421, 6, 28, 29,
                                                                       1806, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5427, 3, 6, 5331,
                                                                       1761, 5337, 1809, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5445, 3, 6, 5337,
                                                                       1764, 5343, 1818, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5463, 3, 6, 5343,
                                                                       1767, 5349, 1827, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5481, 3, 6, 5349,
                                                                       1770, 5355, 1836, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5499, 3, 6, 5355,
                                                                       1773, 5361, 1845, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5517, 3, 6, 5361,
                                                                       1776, 5367, 1854, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5535, 3, 6, 5367,
                                                                       1779, 5373, 1863, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5553, 3, 6, 5379,
                                                                       1785, 5385, 1872, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5571, 3, 6, 5385,
                                                                       1788, 5391, 1881, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5589, 3, 6, 5391,
                                                                       1791, 5397, 1890, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5607, 3, 6, 5397,
                                                                       1794, 5403, 1899, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5625, 3, 6, 5403,
                                                                       1797, 5409, 1908, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5643, 3, 6, 5409,
                                                                       1800, 5415, 1917, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5661, 3, 6, 5415,
                                                                       1803, 5421, 1926, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5679, 0, 6, 5331,
                                                                       1761, 5337, 1935, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5697, 0, 6, 5337,
                                                                       1764, 5343, 1944, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5715, 0, 6, 5343,
                                                                       1767, 5349, 1953, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5733, 0, 6, 5349,
                                                                       1770, 5355, 1962, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5751, 0, 6, 5355,
                                                                       1773, 5361, 1971, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5769, 0, 6, 5361,
                                                                       1776, 5367, 1980, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5787, 0, 6, 5367,
                                                                       1779, 5373, 1989, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5805, 0, 6, 5379,
                                                                       1785, 5385, 1998, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5823, 0, 6, 5385,
                                                                       1788, 5391, 2007, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5841, 0, 6, 5391,
                                                                       1791, 5397, 2016, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5859, 0, 6, 5397,
                                                                       1794, 5403, 2025, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5877, 0, 6, 5403,
                                                                       1797, 5409, 2034, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5895, 0, 6, 5409,
                                                                       1800, 5415, 2043, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5913, 0, 6, 5415,
                                                                       1803, 5421, 2052, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5931, 0, 3, 6,
                                                                       5427, 1809, 5445, 5679,
                                                                       1935, 5697, 139, 148,
                                                                       2061, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5985, 0, 3, 6,
                                                                       5445, 1818, 5463, 5697,
                                                                       1944, 5715, 148, 157,
                                                                       2088, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6039, 0, 3, 6,
                                                                       5463, 1827, 5481, 5715,
                                                                       1953, 5733, 157, 166,
                                                                       2115, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6093, 0, 3, 6,
                                                                       5481, 1836, 5499, 5733,
                                                                       1962, 5751, 166, 175,
                                                                       2142, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6147, 0, 3, 6,
                                                                       5499, 1845, 5517, 5751,
                                                                       1971, 5769, 175, 184,
                                                                       2169, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6201, 0, 3, 6,
                                                                       5517, 1854, 5535, 5769,
                                                                       1980, 5787, 184, 193,
                                                                       2196, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6255, 0, 3, 6,
                                                                       5553, 1872, 5571, 5805,
                                                                       1998, 5823, 211, 220,
                                                                       2223, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6309, 0, 3, 6,
                                                                       5571, 1881, 5589, 5823,
                                                                       2007, 5841, 220, 229,
                                                                       2250, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6363, 0, 3, 6,
                                                                       5589, 1890, 5607, 5841,
                                                                       2016, 5859, 229, 238,
                                                                       2277, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6417, 0, 3, 6,
                                                                       5607, 1899, 5625, 5859,
                                                                       2025, 5877, 238, 247,
                                                                       2304, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6471, 0, 3, 6,
                                                                       5625, 1908, 5643, 5877,
                                                                       2034, 5895, 247, 256,
                                                                       2331, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6525, 0, 3, 6,
                                                                       5643, 1917, 5661, 5895,
                                                                       2043, 5913, 256, 265,
                                                                       2358, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6579, 0, 6, 5679,
                                                                       1935, 5697, 283, 289,
                                                                       2385, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6615, 0, 6, 5697,
                                                                       1944, 5715, 289, 295,
                                                                       2403, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6651, 0, 6, 5715,
                                                                       1953, 5733, 295, 301,
                                                                       2421, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6687, 0, 6, 5733,
                                                                       1962, 5751, 301, 307,
                                                                       2439, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6723, 0, 6, 5751,
                                                                       1971, 5769, 307, 313,
                                                                       2457, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6759, 0, 6, 5769,
                                                                       1980, 5787, 313, 319,
                                                                       2475, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6795, 0, 6, 5805,
                                                                       1998, 5823, 331, 337,
                                                                       2493, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6831, 0, 6, 5823,
                                                                       2007, 5841, 337, 343,
                                                                       2511, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6867, 0, 6, 5841,
                                                                       2016, 5859, 343, 349,
                                                                       2529, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6903, 0, 6, 5859,
                                                                       2025, 5877, 349, 355,
                                                                       2547, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6939, 0, 6, 5877,
                                                                       2034, 5895, 355, 361,
                                                                       2565, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6975, 0, 6, 5895,
                                                                       2043, 5913, 361, 367,
                                                                       2583, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 7011, 0, 3, 6,
                                                                       5931, 2061, 5985, 6579,
                                                                       2385, 6615, 379, 397,
                                                                       2601, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 7119, 0, 3, 6,
                                                                       5985, 2088, 6039, 6615,
                                                                       2403, 6651, 397, 415,
                                                                       2655, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 7227, 0, 3, 6,
                                                                       6039, 2115, 6093, 6651,
                                                                       2421, 6687, 415, 433,
                                                                       2709, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 7335, 0, 3, 6,
                                                                       6093, 2142, 6147, 6687,
                                                                       2439, 6723, 433, 451,
                                                                       2763, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 7443, 0, 3, 6,
                                                                       6147, 2169, 6201, 6723,
                                                                       2457, 6759, 451, 469,
                                                                       2817, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 7551, 0, 3, 6,
                                                                       6255, 2223, 6309, 6795,
                                                                       2493, 6831, 505, 523,
                                                                       2871, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 7659, 0, 3, 6,
                                                                       6309, 2250, 6363, 6831,
                                                                       2511, 6867, 523, 541,
                                                                       2925, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 7767, 0, 3, 6,
                                                                       6363, 2277, 6417, 6867,
                                                                       2529, 6903, 541, 559,
                                                                       2979, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 7875, 0, 3, 6,
                                                                       6417, 2304, 6471, 6903,
                                                                       2547, 6939, 559, 577,
                                                                       3033, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 7983, 0, 3, 6,
                                                                       6471, 2331, 6525, 6939,
                                                                       2565, 6975, 577, 595,
                                                                       3087, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8091, 0, 6, 6579,
                                                                       2385, 6615, 631, 641,
                                                                       3141, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8151, 0, 6, 6615,
                                                                       2403, 6651, 641, 651,
                                                                       3171, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8211, 0, 6, 6651,
                                                                       2421, 6687, 651, 661,
                                                                       3201, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8271, 0, 6, 6687,
                                                                       2439, 6723, 661, 671,
                                                                       3231, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8331, 0, 6, 6723,
                                                                       2457, 6759, 671, 681,
                                                                       3261, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8391, 0, 6, 6795,
                                                                       2493, 6831, 701, 711,
                                                                       3291, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8451, 0, 6, 6831,
                                                                       2511, 6867, 711, 721,
                                                                       3321, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8511, 0, 6, 6867,
                                                                       2529, 6903, 721, 731,
                                                                       3351, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8571, 0, 6, 6903,
                                                                       2547, 6939, 731, 741,
                                                                       3381, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8631, 0, 6, 6939,
                                                                       2565, 6975, 741, 751,
                                                                       3411, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 8691, 0, 3, 6,
                                                                       5931, 5985, 7011, 2601,
                                                                       7119, 8091, 3141, 8151,
                                                                       771, 801, 3441, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 8871, 0, 3, 6,
                                                                       5985, 6039, 7119, 2655,
                                                                       7227, 8151, 3171, 8211,
                                                                       801, 831, 3531, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 9051, 0, 3, 6,
                                                                       6039, 6093, 7227, 2709,
                                                                       7335, 8211, 3201, 8271,
                                                                       831, 861, 3621, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 9231, 0, 3, 6,
                                                                       6093, 6147, 7335, 2763,
                                                                       7443, 8271, 3231, 8331,
                                                                       861, 891, 3711, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 9411, 0, 3, 6,
                                                                       6255, 6309, 7551, 2871,
                                                                       7659, 8391, 3291, 8451,
                                                                       951, 981, 3801, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 9591, 0, 3, 6,
                                                                       6309, 6363, 7659, 2925,
                                                                       7767, 8451, 3321, 8511,
                                                                       981, 1011, 3891, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 9771, 0, 3, 6,
                                                                       6363, 6417, 7767, 2979,
                                                                       7875, 8511, 3351, 8571,
                                                                       1011, 1041, 3981, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 9951, 0, 3, 6,
                                                                       6417, 6471, 7875, 3033,
                                                                       7983, 8571, 3381, 8631,
                                                                       1041, 1071, 4071, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10131, 0, 6, 8091,
                                                                       3141, 8151, 1131, 1146,
                                                                       4161, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10221, 0, 6, 8151,
                                                                       3171, 8211, 1146, 1161,
                                                                       4206, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10311, 0, 6, 8211,
                                                                       3201, 8271, 1161, 1176,
                                                                       4251, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10401, 0, 6, 8271,
                                                                       3231, 8331, 1176, 1191,
                                                                       4296, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10491, 0, 6, 8391,
                                                                       3291, 8451, 1221, 1236,
                                                                       4341, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10581, 0, 6, 8451,
                                                                       3321, 8511, 1236, 1251,
                                                                       4386, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10671, 0, 6, 8511,
                                                                       3351, 8571, 1251, 1266,
                                                                       4431, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 10761, 0, 6, 8571,
                                                                       3381, 8631, 1266, 1281,
                                                                       4476, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 10851, 0, 3, 6,
                                                                       7011, 7119, 8691, 3441,
                                                                       8871, 10131, 4161, 10221,
                                                                       1311, 1356, 4521, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 11121, 0, 3, 6,
                                                                       7119, 7227, 8871, 3531,
                                                                       9051, 10221, 4206, 10311,
                                                                       1356, 1401, 4656, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 11391, 0, 3, 6,
                                                                       7227, 7335, 9051, 3621,
                                                                       9231, 10311, 4251, 10401,
                                                                       1401, 1446, 4791, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 11661, 0, 3, 6,
                                                                       7551, 7659, 9411, 3801,
                                                                       9591, 10491, 4341, 10581,
                                                                       1536, 1581, 4926, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 11931, 0, 3, 6,
                                                                       7659, 7767, 9591, 3891,
                                                                       9771, 10581, 4386, 10671,
                                                                       1581, 1626, 5061, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 12201, 0, 3, 6,
                                                                       7767, 7875, 9771, 3981,
                                                                       9951, 10671, 4431, 10761,
                                                                       1626, 1671, 5196, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12471, 6, 1761,
                                                                       1764, 5343, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12481, 6, 1764,
                                                                       1767, 5349, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12491, 6, 1767,
                                                                       1770, 5355, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12501, 6, 1770,
                                                                       1773, 5361, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12511, 6, 1773,
                                                                       1776, 5367, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12521, 6, 1776,
                                                                       1779, 5373, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12531, 6, 1785,
                                                                       1788, 5391, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12541, 6, 1788,
                                                                       1791, 5397, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12551, 6, 1791,
                                                                       1794, 5403, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12561, 6, 1794,
                                                                       1797, 5409, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12571, 6, 1797,
                                                                       1800, 5415, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12581, 6, 1800,
                                                                       1803, 5421, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12591, 3, 6,
                                                                       12471, 5343, 12481, 5463,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12621, 3, 6,
                                                                       12481, 5349, 12491, 5481,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12651, 3, 6,
                                                                       12491, 5355, 12501, 5499,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12681, 3, 6,
                                                                       12501, 5361, 12511, 5517,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12711, 3, 6,
                                                                       12511, 5367, 12521, 5535,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12741, 3, 6,
                                                                       12531, 5391, 12541, 5589,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12771, 3, 6,
                                                                       12541, 5397, 12551, 5607,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12801, 3, 6,
                                                                       12551, 5403, 12561, 5625,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12831, 3, 6,
                                                                       12561, 5409, 12571, 5643,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12861, 3, 6,
                                                                       12571, 5415, 12581, 5661,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12891, 0, 6,
                                                                       12471, 5343, 12481, 1935,
                                                                       1944, 5715, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12921, 0, 6,
                                                                       12481, 5349, 12491, 1944,
                                                                       1953, 5733, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12951, 0, 6,
                                                                       12491, 5355, 12501, 1953,
                                                                       1962, 5751, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12981, 0, 6,
                                                                       12501, 5361, 12511, 1962,
                                                                       1971, 5769, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13011, 0, 6,
                                                                       12511, 5367, 12521, 1971,
                                                                       1980, 5787, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13041, 0, 6,
                                                                       12531, 5391, 12541, 1998,
                                                                       2007, 5841, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13071, 0, 6,
                                                                       12541, 5397, 12551, 2007,
                                                                       2016, 5859, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13101, 0, 6,
                                                                       12551, 5403, 12561, 2016,
                                                                       2025, 5877, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13131, 0, 6,
                                                                       12561, 5409, 12571, 2025,
                                                                       2034, 5895, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 13161, 0, 6,
                                                                       12571, 5415, 12581, 2034,
                                                                       2043, 5913, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 13191, 0, 3, 6,
                                                                       12591, 5463, 12621, 12891,
                                                                       5715, 12921, 2061, 2088,
                                                                       6039, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 13281, 0, 3, 6,
                                                                       12621, 5481, 12651, 12921,
                                                                       5733, 12951, 2088, 2115,
                                                                       6093, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 13371, 0, 3, 6,
                                                                       12651, 5499, 12681, 12951,
                                                                       5751, 12981, 2115, 2142,
                                                                       6147, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 13461, 0, 3, 6,
                                                                       12681, 5517, 12711, 12981,
                                                                       5769, 13011, 2142, 2169,
                                                                       6201, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 13551, 0, 3, 6,
                                                                       12741, 5589, 12771, 13041,
                                                                       5841, 13071, 2223, 2250,
                                                                       6363, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 13641, 0, 3, 6,
                                                                       12771, 5607, 12801, 13071,
                                                                       5859, 13101, 2250, 2277,
                                                                       6417, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 13731, 0, 3, 6,
                                                                       12801, 5625, 12831, 13101,
                                                                       5877, 13131, 2277, 2304,
                                                                       6471, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 13821, 0, 3, 6,
                                                                       12831, 5643, 12861, 13131,
                                                                       5895, 13161, 2304, 2331,
                                                                       6525, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13911, 0, 6,
                                                                       12891, 5715, 12921, 2385,
                                                                       2403, 6651, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 13971, 0, 6,
                                                                       12921, 5733, 12951, 2403,
                                                                       2421, 6687, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14031, 0, 6,
                                                                       12951, 5751, 12981, 2421,
                                                                       2439, 6723, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14091, 0, 6,
                                                                       12981, 5769, 13011, 2439,
                                                                       2457, 6759, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14151, 0, 6,
                                                                       13041, 5841, 13071, 2493,
                                                                       2511, 6867, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14211, 0, 6,
                                                                       13071, 5859, 13101, 2511,
                                                                       2529, 6903, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14271, 0, 6,
                                                                       13101, 5877, 13131, 2529,
                                                                       2547, 6939, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 14331, 0, 6,
                                                                       13131, 5895, 13161, 2547,
                                                                       2565, 6975, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 14391, 0, 3, 6,
                                                                       13191, 6039, 13281, 13911,
                                                                       6651, 13971, 2601, 2655,
                                                                       7227, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 14571, 0, 3, 6,
                                                                       13281, 6093, 13371, 13971,
                                                                       6687, 14031, 2655, 2709,
                                                                       7335, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 14751, 0, 3, 6,
                                                                       13371, 6147, 13461, 14031,
                                                                       6723, 14091, 2709, 2763,
                                                                       7443, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 14931, 0, 3, 6,
                                                                       13551, 6363, 13641, 14151,
                                                                       6867, 14211, 2871, 2925,
                                                                       7767, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 15111, 0, 3, 6,
                                                                       13641, 6417, 13731, 14211,
                                                                       6903, 14271, 2925, 2979,
                                                                       7875, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 15291, 0, 3, 6,
                                                                       13731, 6471, 13821, 14271,
                                                                       6939, 14331, 2979, 3033,
                                                                       7983, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15471, 0, 6,
                                                                       13911, 6651, 13971, 3141,
                                                                       3171, 8211, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15571, 0, 6,
                                                                       13971, 6687, 14031, 3171,
                                                                       3201, 8271, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15671, 0, 6,
                                                                       14031, 6723, 14091, 3201,
                                                                       3231, 8331, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15771, 0, 6,
                                                                       14151, 6867, 14211, 3291,
                                                                       3321, 8511, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15871, 0, 6,
                                                                       14211, 6903, 14271, 3321,
                                                                       3351, 8571, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 15971, 0, 6,
                                                                       14271, 6939, 14331, 3351,
                                                                       3381, 8631, ncols, gamma,
                                                                       p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 16071, 0, 3, 6,
                                                                       13191, 13281, 14391, 7227,
                                                                       14571, 15471, 8211, 15571,
                                                                       3441, 3531, 9051, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 16371, 0, 3, 6,
                                                                       13281, 13371, 14571, 7335,
                                                                       14751, 15571, 8271, 15671,
                                                                       3531, 3621, 9231, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 16671, 0, 3, 6,
                                                                       13551, 13641, 14931, 7767,
                                                                       15111, 15771, 8511, 15871,
                                                                       3801, 3891, 9771, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 16971, 0, 3, 6,
                                                                       13641, 13731, 15111, 7875,
                                                                       15291, 15871, 8571, 15971,
                                                                       3891, 3981, 9951, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17271, 0, 6,
                                                                       15471, 8211, 15571, 4161,
                                                                       4206, 10311, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17421, 0, 6,
                                                                       15571, 8271, 15671, 4206,
                                                                       4251, 10401, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17571, 0, 6,
                                                                       15771, 8511, 15871, 4341,
                                                                       4386, 10671, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 17721, 0, 6,
                                                                       15871, 8571, 15971, 4386,
                                                                       4431, 10761, ncols, gamma,
                                                                       p, q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 17871, 0, 3, 6,
                                                                       14391, 14571, 16071, 9051,
                                                                       16371, 17271, 10311,
                                                                       17421, 4521, 4656, 11391,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 18321, 0, 3, 6,
                                                                       14931, 15111, 16671, 9771,
                                                                       16971, 17571, 10671,
                                                                       17721, 4926, 5061, 12201,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18771, 6, 5331,
                                                                       5337, 12471, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18786, 6, 5337,
                                                                       5343, 12481, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18801, 6, 5343,
                                                                       5349, 12491, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18816, 6, 5349,
                                                                       5355, 12501, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18831, 6, 5355,
                                                                       5361, 12511, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18846, 6, 5361,
                                                                       5367, 12521, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18861, 6, 5379,
                                                                       5385, 12531, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18876, 6, 5385,
                                                                       5391, 12541, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18891, 6, 5391,
                                                                       5397, 12551, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18906, 6, 5397,
                                                                       5403, 12561, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18921, 6, 5403,
                                                                       5409, 12571, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18936, 6, 5409,
                                                                       5415, 12581, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18951, 3, 6,
                                                                       18771, 12471, 18786, 5427,
                                                                       5445, 12591, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18996, 3, 6,
                                                                       18786, 12481, 18801, 5445,
                                                                       5463, 12621, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19041, 3, 6,
                                                                       18801, 12491, 18816, 5463,
                                                                       5481, 12651, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19086, 3, 6,
                                                                       18816, 12501, 18831, 5481,
                                                                       5499, 12681, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19131, 3, 6,
                                                                       18831, 12511, 18846, 5499,
                                                                       5517, 12711, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19176, 3, 6,
                                                                       18861, 12531, 18876, 5553,
                                                                       5571, 12741, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19221, 3, 6,
                                                                       18876, 12541, 18891, 5571,
                                                                       5589, 12771, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19266, 3, 6,
                                                                       18891, 12551, 18906, 5589,
                                                                       5607, 12801, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19311, 3, 6,
                                                                       18906, 12561, 18921, 5607,
                                                                       5625, 12831, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19356, 3, 6,
                                                                       18921, 12571, 18936, 5625,
                                                                       5643, 12861, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19401, 0, 6,
                                                                       18771, 12471, 18786, 5679,
                                                                       5697, 12891, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19446, 0, 6,
                                                                       18786, 12481, 18801, 5697,
                                                                       5715, 12921, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19491, 0, 6,
                                                                       18801, 12491, 18816, 5715,
                                                                       5733, 12951, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19536, 0, 6,
                                                                       18816, 12501, 18831, 5733,
                                                                       5751, 12981, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19581, 0, 6,
                                                                       18831, 12511, 18846, 5751,
                                                                       5769, 13011, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19626, 0, 6,
                                                                       18861, 12531, 18876, 5805,
                                                                       5823, 13041, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19671, 0, 6,
                                                                       18876, 12541, 18891, 5823,
                                                                       5841, 13071, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19716, 0, 6,
                                                                       18891, 12551, 18906, 5841,
                                                                       5859, 13101, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19761, 0, 6,
                                                                       18906, 12561, 18921, 5859,
                                                                       5877, 13131, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19806, 0, 6,
                                                                       18921, 12571, 18936, 5877,
                                                                       5895, 13161, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 19851, 0, 3, 6,
                                                                       18951, 12591, 18996,
                                                                       19401, 12891, 19446, 5931,
                                                                       5985, 13191, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 19986, 0, 3, 6,
                                                                       18996, 12621, 19041,
                                                                       19446, 12921, 19491, 5985,
                                                                       6039, 13281, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 20121, 0, 3, 6,
                                                                       19041, 12651, 19086,
                                                                       19491, 12951, 19536, 6039,
                                                                       6093, 13371, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 20256, 0, 3, 6,
                                                                       19086, 12681, 19131,
                                                                       19536, 12981, 19581, 6093,
                                                                       6147, 13461, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 20391, 0, 3, 6,
                                                                       19176, 12741, 19221,
                                                                       19626, 13041, 19671, 6255,
                                                                       6309, 13551, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 20526, 0, 3, 6,
                                                                       19221, 12771, 19266,
                                                                       19671, 13071, 19716, 6309,
                                                                       6363, 13641, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 20661, 0, 3, 6,
                                                                       19266, 12801, 19311,
                                                                       19716, 13101, 19761, 6363,
                                                                       6417, 13731, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 20796, 0, 3, 6,
                                                                       19311, 12831, 19356,
                                                                       19761, 13131, 19806, 6417,
                                                                       6471, 13821, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 20931, 0, 6,
                                                                       19401, 12891, 19446, 6579,
                                                                       6615, 13911, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21021, 0, 6,
                                                                       19446, 12921, 19491, 6615,
                                                                       6651, 13971, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21111, 0, 6,
                                                                       19491, 12951, 19536, 6651,
                                                                       6687, 14031, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21201, 0, 6,
                                                                       19536, 12981, 19581, 6687,
                                                                       6723, 14091, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21291, 0, 6,
                                                                       19626, 13041, 19671, 6795,
                                                                       6831, 14151, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21381, 0, 6,
                                                                       19671, 13071, 19716, 6831,
                                                                       6867, 14211, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21471, 0, 6,
                                                                       19716, 13101, 19761, 6867,
                                                                       6903, 14271, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 21561, 0, 6,
                                                                       19761, 13131, 19806, 6903,
                                                                       6939, 14331, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 21651, 0, 3, 6,
                                                                       19851, 13191, 19986,
                                                                       20931, 13911, 21021, 7011,
                                                                       7119, 14391, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 21921, 0, 3, 6,
                                                                       19986, 13281, 20121,
                                                                       21021, 13971, 21111, 7119,
                                                                       7227, 14571, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 22191, 0, 3, 6,
                                                                       20121, 13371, 20256,
                                                                       21111, 14031, 21201, 7227,
                                                                       7335, 14751, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 22461, 0, 3, 6,
                                                                       20391, 13551, 20526,
                                                                       21291, 14151, 21381, 7551,
                                                                       7659, 14931, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 22731, 0, 3, 6,
                                                                       20526, 13641, 20661,
                                                                       21381, 14211, 21471, 7659,
                                                                       7767, 15111, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 23001, 0, 3, 6,
                                                                       20661, 13731, 20796,
                                                                       21471, 14271, 21561, 7767,
                                                                       7875, 15291, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 23271, 0, 6,
                                                                       20931, 13911, 21021, 8091,
                                                                       8151, 15471, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 23421, 0, 6,
                                                                       21021, 13971, 21111, 8151,
                                                                       8211, 15571, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 23571, 0, 6,
                                                                       21111, 14031, 21201, 8211,
                                                                       8271, 15671, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 23721, 0, 6,
                                                                       21291, 14151, 21381, 8391,
                                                                       8451, 15771, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 23871, 0, 6,
                                                                       21381, 14211, 21471, 8451,
                                                                       8511, 15871, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 24021, 0, 6,
                                                                       21471, 14271, 21561, 8511,
                                                                       8571, 15971, ncols, gamma,
                                                                       p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 24171, 0, 3, 6,
                                                                       19851, 19986, 21651,
                                                                       14391, 21921, 23271,
                                                                       15471, 23421, 8691, 8871,
                                                                       16071, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 24621, 0, 3, 6,
                                                                       19986, 20121, 21921,
                                                                       14571, 22191, 23421,
                                                                       15571, 23571, 8871, 9051,
                                                                       16371, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 25071, 0, 3, 6,
                                                                       20391, 20526, 22461,
                                                                       14931, 22731, 23721,
                                                                       15771, 23871, 9411, 9591,
                                                                       16671, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 25521, 0, 3, 6,
                                                                       20526, 20661, 22731,
                                                                       15111, 23001, 23871,
                                                                       15871, 24021, 9591, 9771,
                                                                       16971, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 25971, 0, 6,
                                                                       23271, 15471, 23421,
                                                                       10131, 10221, 17271,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 26196, 0, 6,
                                                                       23421, 15571, 23571,
                                                                       10221, 10311, 17421,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 26421, 0, 6,
                                                                       23721, 15771, 23871,
                                                                       10491, 10581, 17571,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 26646, 0, 6,
                                                                       23871, 15871, 24021,
                                                                       10581, 10671, 17721,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 26871, 0, 3, 6,
                                                                       21651, 21921, 24171,
                                                                       16071, 24621, 25971,
                                                                       17271, 26196, 10851,
                                                                       11121, 17871, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 27546, 0, 3, 6,
                                                                       22461, 22731, 25071,
                                                                       16671, 25521, 26421,
                                                                       17571, 26646, 11661,
                                                                       11931, 18321, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 28221, 27546, 15, 15, ncols, beta);

                    simdgeo::geom_s_y(buffer, 28446, 27546, 15, 15, ncols, beta);

                    simdgeo::geom_s_z(buffer, 28671, 27546, 15, 15, ncols, beta);

                    simdgeo::geom_s_x(buffer, 28896, 26871, 15, 15, ncols, beta);

                    simdgeo::geom_s_y(buffer, 29121, 26871, 15, 15, ncols, beta);

                    simdgeo::geom_s_z(buffer, 29346, 26871, 15, 15, ncols, beta);

                    simdfunc::contract_primitives(buffer, 29571, 28221, 1350, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 30921, 29571, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 30921, 9, nmax);

        simdtrf::transform_g_inner(buffer, 30921, 29796, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 81 * nvalues + n * npairs, nvalues, buffer, 30921, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 30921, 30021, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 162 * nvalues + n * npairs, nvalues, buffer, 30921,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 30921, 30246, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 243 * nvalues + n * npairs, nvalues, buffer, 30921,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 30921, 30471, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 324 * nvalues + n * npairs, nvalues, buffer, 30921,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 30921, 30696, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 405 * nvalues + n * npairs, nvalues, buffer, 30921,
                                   9, nmax);
    }

    for (size_t m = 0; m < 486; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
