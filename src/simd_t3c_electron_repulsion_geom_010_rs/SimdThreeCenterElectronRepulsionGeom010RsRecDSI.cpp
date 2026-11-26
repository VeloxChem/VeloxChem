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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecDSI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_dsi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_dsi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 21225, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 390 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 21225, 20139, 1008, dimensions);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 631, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 634, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 637, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 640, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 643, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 646, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 649, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 652, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 655, 6, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 658, 6, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 661, 6, 25, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 664, 6, 26, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 667, 6, 27, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 670, 6, 28, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 673, 6, 29, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 676, 6, 30, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 679, 6, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 688, 6, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 697, 6, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 706, 6, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 715, 6, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 724, 6, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 733, 6, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 742, 6, 23, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 751, 6, 24, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 760, 6, 25, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 769, 6, 26, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 778, 6, 27, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 787, 6, 28, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 796, 6, 29, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 805, 6, 12, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 814, 6, 13, 94,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 823, 6, 14, 97,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 832, 6, 15, 100,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 841, 6, 16, 103,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 850, 6, 17, 106,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 859, 6, 18, 109,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 868, 6, 23, 118,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 877, 6, 24, 121,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 886, 6, 25, 124,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 895, 6, 26, 127,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 904, 6, 27, 130,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 913, 6, 28, 133,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 922, 6, 29, 136,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 931, 6, 37, 91,
                                                                       157, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 958, 6, 40, 94,
                                                                       166, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 985, 6, 43, 97,
                                                                       175, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1012, 6, 46, 100,
                                                                       184, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1039, 6, 49, 103,
                                                                       193, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1066, 6, 52, 106,
                                                                       202, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1093, 6, 64, 118,
                                                                       229, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1120, 6, 67, 121,
                                                                       238, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1147, 6, 70, 124,
                                                                       247, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1174, 6, 73, 127,
                                                                       256, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1201, 6, 76, 130,
                                                                       265, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1228, 6, 79, 133,
                                                                       274, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1255, 6, 91, 295,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1273, 6, 94, 301,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1291, 6, 97, 307,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1309, 6, 100, 313,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1327, 6, 103, 319,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1345, 6, 106, 325,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1363, 6, 118, 343,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1381, 6, 121, 349,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1399, 6, 124, 355,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1417, 6, 127, 361,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1435, 6, 130, 367,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1453, 6, 133, 373,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1471, 0, 6, 931,
                                                                       157, 958, 295, 415, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1525, 0, 6, 958,
                                                                       166, 985, 301, 433, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1579, 0, 6, 985,
                                                                       175, 1012, 307, 451,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1633, 0, 6, 1012,
                                                                       184, 1039, 313, 469,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1687, 0, 6, 1039,
                                                                       193, 1066, 319, 487,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1741, 0, 6, 1093,
                                                                       229, 1120, 343, 541,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1795, 0, 6, 1120,
                                                                       238, 1147, 349, 559,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1849, 0, 6, 1147,
                                                                       247, 1174, 355, 577,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1903, 0, 6, 1174,
                                                                       256, 1201, 361, 595,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1957, 0, 6, 1201,
                                                                       265, 1228, 367, 613,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2011, 6, 10, 11,
                                                                       631, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2017, 6, 11, 12,
                                                                       634, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2023, 6, 12, 13,
                                                                       637, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2029, 6, 13, 14,
                                                                       640, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2035, 6, 14, 15,
                                                                       643, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2041, 6, 15, 16,
                                                                       646, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2047, 6, 16, 17,
                                                                       649, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2053, 6, 17, 18,
                                                                       652, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2059, 6, 21, 22,
                                                                       655, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2065, 6, 22, 23,
                                                                       658, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2071, 6, 23, 24,
                                                                       661, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2077, 6, 24, 25,
                                                                       664, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2083, 6, 25, 26,
                                                                       667, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2089, 6, 26, 27,
                                                                       670, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2095, 6, 27, 28,
                                                                       673, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2101, 6, 28, 29,
                                                                       676, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2107, 3, 6, 2011,
                                                                       631, 2017, 679, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2125, 3, 6, 2017,
                                                                       634, 2023, 688, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2143, 3, 6, 2023,
                                                                       637, 2029, 697, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2161, 3, 6, 2029,
                                                                       640, 2035, 706, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2179, 3, 6, 2035,
                                                                       643, 2041, 715, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2197, 3, 6, 2041,
                                                                       646, 2047, 724, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2215, 3, 6, 2047,
                                                                       649, 2053, 733, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2233, 3, 6, 2059,
                                                                       655, 2065, 742, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2251, 3, 6, 2065,
                                                                       658, 2071, 751, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2269, 3, 6, 2071,
                                                                       661, 2077, 760, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2287, 3, 6, 2077,
                                                                       664, 2083, 769, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2305, 3, 6, 2083,
                                                                       667, 2089, 778, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2323, 3, 6, 2089,
                                                                       670, 2095, 787, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2341, 3, 6, 2095,
                                                                       673, 2101, 796, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2359, 0, 6, 2011,
                                                                       631, 2017, 805, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2377, 0, 6, 2017,
                                                                       634, 2023, 814, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2395, 0, 6, 2023,
                                                                       637, 2029, 823, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2413, 0, 6, 2029,
                                                                       640, 2035, 832, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2431, 0, 6, 2035,
                                                                       643, 2041, 841, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2449, 0, 6, 2041,
                                                                       646, 2047, 850, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2467, 0, 6, 2047,
                                                                       649, 2053, 859, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2485, 0, 6, 2059,
                                                                       655, 2065, 868, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2503, 0, 6, 2065,
                                                                       658, 2071, 877, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2521, 0, 6, 2071,
                                                                       661, 2077, 886, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2539, 0, 6, 2077,
                                                                       664, 2083, 895, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2557, 0, 6, 2083,
                                                                       667, 2089, 904, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2575, 0, 6, 2089,
                                                                       670, 2095, 913, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2593, 0, 6, 2095,
                                                                       673, 2101, 922, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2611, 0, 3, 6,
                                                                       2107, 679, 2125, 2359,
                                                                       805, 2377, 139, 148, 931,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2665, 0, 3, 6,
                                                                       2125, 688, 2143, 2377,
                                                                       814, 2395, 148, 157, 958,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2719, 0, 3, 6,
                                                                       2143, 697, 2161, 2395,
                                                                       823, 2413, 157, 166, 985,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2773, 0, 3, 6,
                                                                       2161, 706, 2179, 2413,
                                                                       832, 2431, 166, 175, 1012,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2827, 0, 3, 6,
                                                                       2179, 715, 2197, 2431,
                                                                       841, 2449, 175, 184, 1039,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2881, 0, 3, 6,
                                                                       2197, 724, 2215, 2449,
                                                                       850, 2467, 184, 193, 1066,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2935, 0, 3, 6,
                                                                       2233, 742, 2251, 2485,
                                                                       868, 2503, 211, 220, 1093,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2989, 0, 3, 6,
                                                                       2251, 751, 2269, 2503,
                                                                       877, 2521, 220, 229, 1120,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3043, 0, 3, 6,
                                                                       2269, 760, 2287, 2521,
                                                                       886, 2539, 229, 238, 1147,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3097, 0, 3, 6,
                                                                       2287, 769, 2305, 2539,
                                                                       895, 2557, 238, 247, 1174,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3151, 0, 3, 6,
                                                                       2305, 778, 2323, 2557,
                                                                       904, 2575, 247, 256, 1201,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3205, 0, 3, 6,
                                                                       2323, 787, 2341, 2575,
                                                                       913, 2593, 256, 265, 1228,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3259, 0, 6, 2359,
                                                                       805, 2377, 283, 289, 1255,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3295, 0, 6, 2377,
                                                                       814, 2395, 289, 295, 1273,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3331, 0, 6, 2395,
                                                                       823, 2413, 295, 301, 1291,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3367, 0, 6, 2413,
                                                                       832, 2431, 301, 307, 1309,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3403, 0, 6, 2431,
                                                                       841, 2449, 307, 313, 1327,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3439, 0, 6, 2449,
                                                                       850, 2467, 313, 319, 1345,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3475, 0, 6, 2485,
                                                                       868, 2503, 331, 337, 1363,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3511, 0, 6, 2503,
                                                                       877, 2521, 337, 343, 1381,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3547, 0, 6, 2521,
                                                                       886, 2539, 343, 349, 1399,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3583, 0, 6, 2539,
                                                                       895, 2557, 349, 355, 1417,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3619, 0, 6, 2557,
                                                                       904, 2575, 355, 361, 1435,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3655, 0, 6, 2575,
                                                                       913, 2593, 361, 367, 1453,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3691, 0, 3, 6,
                                                                       2611, 931, 2665, 3259,
                                                                       1255, 3295, 379, 397,
                                                                       1471, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3799, 0, 3, 6,
                                                                       2665, 958, 2719, 3295,
                                                                       1273, 3331, 397, 415,
                                                                       1525, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3907, 0, 3, 6,
                                                                       2719, 985, 2773, 3331,
                                                                       1291, 3367, 415, 433,
                                                                       1579, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4015, 0, 3, 6,
                                                                       2773, 1012, 2827, 3367,
                                                                       1309, 3403, 433, 451,
                                                                       1633, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4123, 0, 3, 6,
                                                                       2827, 1039, 2881, 3403,
                                                                       1327, 3439, 451, 469,
                                                                       1687, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4231, 0, 3, 6,
                                                                       2935, 1093, 2989, 3475,
                                                                       1363, 3511, 505, 523,
                                                                       1741, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4339, 0, 3, 6,
                                                                       2989, 1120, 3043, 3511,
                                                                       1381, 3547, 523, 541,
                                                                       1795, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4447, 0, 3, 6,
                                                                       3043, 1147, 3097, 3547,
                                                                       1399, 3583, 541, 559,
                                                                       1849, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4555, 0, 3, 6,
                                                                       3097, 1174, 3151, 3583,
                                                                       1417, 3619, 559, 577,
                                                                       1903, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4663, 0, 3, 6,
                                                                       3151, 1201, 3205, 3619,
                                                                       1435, 3655, 577, 595,
                                                                       1957, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4771, 6, 631, 634,
                                                                       2023, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4781, 6, 634, 637,
                                                                       2029, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4791, 6, 637, 640,
                                                                       2035, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4801, 6, 640, 643,
                                                                       2041, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4811, 6, 643, 646,
                                                                       2047, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4821, 6, 646, 649,
                                                                       2053, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4831, 6, 655, 658,
                                                                       2071, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4841, 6, 658, 661,
                                                                       2077, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4851, 6, 661, 664,
                                                                       2083, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4861, 6, 664, 667,
                                                                       2089, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4871, 6, 667, 670,
                                                                       2095, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4881, 6, 670, 673,
                                                                       2101, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4891, 3, 6, 4771,
                                                                       2023, 4781, 2143, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4921, 3, 6, 4781,
                                                                       2029, 4791, 2161, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4951, 3, 6, 4791,
                                                                       2035, 4801, 2179, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4981, 3, 6, 4801,
                                                                       2041, 4811, 2197, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5011, 3, 6, 4811,
                                                                       2047, 4821, 2215, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5041, 3, 6, 4831,
                                                                       2071, 4841, 2269, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5071, 3, 6, 4841,
                                                                       2077, 4851, 2287, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5101, 3, 6, 4851,
                                                                       2083, 4861, 2305, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5131, 3, 6, 4861,
                                                                       2089, 4871, 2323, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5161, 3, 6, 4871,
                                                                       2095, 4881, 2341, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5191, 0, 6, 4771,
                                                                       2023, 4781, 805, 814,
                                                                       2395, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5221, 0, 6, 4781,
                                                                       2029, 4791, 814, 823,
                                                                       2413, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5251, 0, 6, 4791,
                                                                       2035, 4801, 823, 832,
                                                                       2431, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5281, 0, 6, 4801,
                                                                       2041, 4811, 832, 841,
                                                                       2449, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5311, 0, 6, 4811,
                                                                       2047, 4821, 841, 850,
                                                                       2467, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5341, 0, 6, 4831,
                                                                       2071, 4841, 868, 877,
                                                                       2521, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5371, 0, 6, 4841,
                                                                       2077, 4851, 877, 886,
                                                                       2539, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5401, 0, 6, 4851,
                                                                       2083, 4861, 886, 895,
                                                                       2557, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5431, 0, 6, 4861,
                                                                       2089, 4871, 895, 904,
                                                                       2575, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5461, 0, 6, 4871,
                                                                       2095, 4881, 904, 913,
                                                                       2593, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5491, 0, 3, 6,
                                                                       4891, 2143, 4921, 5191,
                                                                       2395, 5221, 931, 958,
                                                                       2719, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5581, 0, 3, 6,
                                                                       4921, 2161, 4951, 5221,
                                                                       2413, 5251, 958, 985,
                                                                       2773, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5671, 0, 3, 6,
                                                                       4951, 2179, 4981, 5251,
                                                                       2431, 5281, 985, 1012,
                                                                       2827, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5761, 0, 3, 6,
                                                                       4981, 2197, 5011, 5281,
                                                                       2449, 5311, 1012, 1039,
                                                                       2881, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5851, 0, 3, 6,
                                                                       5041, 2269, 5071, 5341,
                                                                       2521, 5371, 1093, 1120,
                                                                       3043, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5941, 0, 3, 6,
                                                                       5071, 2287, 5101, 5371,
                                                                       2539, 5401, 1120, 1147,
                                                                       3097, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6031, 0, 3, 6,
                                                                       5101, 2305, 5131, 5401,
                                                                       2557, 5431, 1147, 1174,
                                                                       3151, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6121, 0, 3, 6,
                                                                       5131, 2323, 5161, 5431,
                                                                       2575, 5461, 1174, 1201,
                                                                       3205, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6211, 0, 6, 5191,
                                                                       2395, 5221, 1255, 1273,
                                                                       3331, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6271, 0, 6, 5221,
                                                                       2413, 5251, 1273, 1291,
                                                                       3367, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6331, 0, 6, 5251,
                                                                       2431, 5281, 1291, 1309,
                                                                       3403, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6391, 0, 6, 5281,
                                                                       2449, 5311, 1309, 1327,
                                                                       3439, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6451, 0, 6, 5341,
                                                                       2521, 5371, 1363, 1381,
                                                                       3547, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6511, 0, 6, 5371,
                                                                       2539, 5401, 1381, 1399,
                                                                       3583, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6571, 0, 6, 5401,
                                                                       2557, 5431, 1399, 1417,
                                                                       3619, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 6631, 0, 6, 5431,
                                                                       2575, 5461, 1417, 1435,
                                                                       3655, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 6691, 0, 3, 6,
                                                                       5491, 2719, 5581, 6211,
                                                                       3331, 6271, 1471, 1525,
                                                                       3907, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 6871, 0, 3, 6,
                                                                       5581, 2773, 5671, 6271,
                                                                       3367, 6331, 1525, 1579,
                                                                       4015, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 7051, 0, 3, 6,
                                                                       5671, 2827, 5761, 6331,
                                                                       3403, 6391, 1579, 1633,
                                                                       4123, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 7231, 0, 3, 6,
                                                                       5851, 3043, 5941, 6451,
                                                                       3547, 6511, 1741, 1795,
                                                                       4447, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 7411, 0, 3, 6,
                                                                       5941, 3097, 6031, 6511,
                                                                       3583, 6571, 1795, 1849,
                                                                       4555, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 7591, 0, 3, 6,
                                                                       6031, 3151, 6121, 6571,
                                                                       3619, 6631, 1849, 1903,
                                                                       4663, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7771, 6, 2011,
                                                                       2017, 4771, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7786, 6, 2017,
                                                                       2023, 4781, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7801, 6, 2023,
                                                                       2029, 4791, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7816, 6, 2029,
                                                                       2035, 4801, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7831, 6, 2035,
                                                                       2041, 4811, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7846, 6, 2041,
                                                                       2047, 4821, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7861, 6, 2059,
                                                                       2065, 4831, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7876, 6, 2065,
                                                                       2071, 4841, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7891, 6, 2071,
                                                                       2077, 4851, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7906, 6, 2077,
                                                                       2083, 4861, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7921, 6, 2083,
                                                                       2089, 4871, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7936, 6, 2089,
                                                                       2095, 4881, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7951, 3, 6, 7771,
                                                                       4771, 7786, 2107, 2125,
                                                                       4891, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7996, 3, 6, 7786,
                                                                       4781, 7801, 2125, 2143,
                                                                       4921, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8041, 3, 6, 7801,
                                                                       4791, 7816, 2143, 2161,
                                                                       4951, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8086, 3, 6, 7816,
                                                                       4801, 7831, 2161, 2179,
                                                                       4981, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8131, 3, 6, 7831,
                                                                       4811, 7846, 2179, 2197,
                                                                       5011, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8176, 3, 6, 7861,
                                                                       4831, 7876, 2233, 2251,
                                                                       5041, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8221, 3, 6, 7876,
                                                                       4841, 7891, 2251, 2269,
                                                                       5071, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8266, 3, 6, 7891,
                                                                       4851, 7906, 2269, 2287,
                                                                       5101, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8311, 3, 6, 7906,
                                                                       4861, 7921, 2287, 2305,
                                                                       5131, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8356, 3, 6, 7921,
                                                                       4871, 7936, 2305, 2323,
                                                                       5161, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8401, 0, 6, 7771,
                                                                       4771, 7786, 2359, 2377,
                                                                       5191, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8446, 0, 6, 7786,
                                                                       4781, 7801, 2377, 2395,
                                                                       5221, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8491, 0, 6, 7801,
                                                                       4791, 7816, 2395, 2413,
                                                                       5251, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8536, 0, 6, 7816,
                                                                       4801, 7831, 2413, 2431,
                                                                       5281, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8581, 0, 6, 7831,
                                                                       4811, 7846, 2431, 2449,
                                                                       5311, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8626, 0, 6, 7861,
                                                                       4831, 7876, 2485, 2503,
                                                                       5341, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8671, 0, 6, 7876,
                                                                       4841, 7891, 2503, 2521,
                                                                       5371, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8716, 0, 6, 7891,
                                                                       4851, 7906, 2521, 2539,
                                                                       5401, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8761, 0, 6, 7906,
                                                                       4861, 7921, 2539, 2557,
                                                                       5431, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 8806, 0, 6, 7921,
                                                                       4871, 7936, 2557, 2575,
                                                                       5461, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 8851, 0, 3, 6,
                                                                       7951, 4891, 7996, 8401,
                                                                       5191, 8446, 2611, 2665,
                                                                       5491, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 8986, 0, 3, 6,
                                                                       7996, 4921, 8041, 8446,
                                                                       5221, 8491, 2665, 2719,
                                                                       5581, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9121, 0, 3, 6,
                                                                       8041, 4951, 8086, 8491,
                                                                       5251, 8536, 2719, 2773,
                                                                       5671, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9256, 0, 3, 6,
                                                                       8086, 4981, 8131, 8536,
                                                                       5281, 8581, 2773, 2827,
                                                                       5761, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9391, 0, 3, 6,
                                                                       8176, 5041, 8221, 8626,
                                                                       5341, 8671, 2935, 2989,
                                                                       5851, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9526, 0, 3, 6,
                                                                       8221, 5071, 8266, 8671,
                                                                       5371, 8716, 2989, 3043,
                                                                       5941, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9661, 0, 3, 6,
                                                                       8266, 5101, 8311, 8716,
                                                                       5401, 8761, 3043, 3097,
                                                                       6031, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9796, 0, 3, 6,
                                                                       8311, 5131, 8356, 8761,
                                                                       5431, 8806, 3097, 3151,
                                                                       6121, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 9931, 0, 6, 8401,
                                                                       5191, 8446, 3259, 3295,
                                                                       6211, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10021, 0, 6, 8446,
                                                                       5221, 8491, 3295, 3331,
                                                                       6271, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10111, 0, 6, 8491,
                                                                       5251, 8536, 3331, 3367,
                                                                       6331, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10201, 0, 6, 8536,
                                                                       5281, 8581, 3367, 3403,
                                                                       6391, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10291, 0, 6, 8626,
                                                                       5341, 8671, 3475, 3511,
                                                                       6451, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10381, 0, 6, 8671,
                                                                       5371, 8716, 3511, 3547,
                                                                       6511, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10471, 0, 6, 8716,
                                                                       5401, 8761, 3547, 3583,
                                                                       6571, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 10561, 0, 6, 8761,
                                                                       5431, 8806, 3583, 3619,
                                                                       6631, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 10651, 0, 3, 6,
                                                                       8851, 5491, 8986, 9931,
                                                                       6211, 10021, 3691, 3799,
                                                                       6691, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 10921, 0, 3, 6,
                                                                       8986, 5581, 9121, 10021,
                                                                       6271, 10111, 3799, 3907,
                                                                       6871, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 11191, 0, 3, 6,
                                                                       9121, 5671, 9256, 10111,
                                                                       6331, 10201, 3907, 4015,
                                                                       7051, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 11461, 0, 3, 6,
                                                                       9391, 5851, 9526, 10291,
                                                                       6451, 10381, 4231, 4339,
                                                                       7231, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 11731, 0, 3, 6,
                                                                       9526, 5941, 9661, 10381,
                                                                       6511, 10471, 4339, 4447,
                                                                       7411, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 12001, 0, 3, 6,
                                                                       9661, 6031, 9796, 10471,
                                                                       6571, 10561, 4447, 4555,
                                                                       7591, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12271, 6, 4771,
                                                                       4781, 7801, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12292, 6, 4781,
                                                                       4791, 7816, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12313, 6, 4791,
                                                                       4801, 7831, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12334, 6, 4801,
                                                                       4811, 7846, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12355, 6, 4831,
                                                                       4841, 7891, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12376, 6, 4841,
                                                                       4851, 7906, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12397, 6, 4851,
                                                                       4861, 7921, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12418, 6, 4861,
                                                                       4871, 7936, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 12439, 3, 6,
                                                                       12271, 7801, 12292, 4891,
                                                                       4921, 8041, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 12502, 3, 6,
                                                                       12292, 7816, 12313, 4921,
                                                                       4951, 8086, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 12565, 3, 6,
                                                                       12313, 7831, 12334, 4951,
                                                                       4981, 8131, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 12628, 3, 6,
                                                                       12355, 7891, 12376, 5041,
                                                                       5071, 8266, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 12691, 3, 6,
                                                                       12376, 7906, 12397, 5071,
                                                                       5101, 8311, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 12754, 3, 6,
                                                                       12397, 7921, 12418, 5101,
                                                                       5131, 8356, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 12817, 0, 6,
                                                                       12271, 7801, 12292, 5191,
                                                                       5221, 8491, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 12880, 0, 6,
                                                                       12292, 7816, 12313, 5221,
                                                                       5251, 8536, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 12943, 0, 6,
                                                                       12313, 7831, 12334, 5251,
                                                                       5281, 8581, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13006, 0, 6,
                                                                       12355, 7891, 12376, 5341,
                                                                       5371, 8716, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13069, 0, 6,
                                                                       12376, 7906, 12397, 5371,
                                                                       5401, 8761, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13132, 0, 6,
                                                                       12397, 7921, 12418, 5401,
                                                                       5431, 8806, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 13195, 0, 3, 6,
                                                                       12439, 8041, 12502, 12817,
                                                                       8491, 12880, 5491, 5581,
                                                                       9121, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 13384, 0, 3, 6,
                                                                       12502, 8086, 12565, 12880,
                                                                       8536, 12943, 5581, 5671,
                                                                       9256, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 13573, 0, 3, 6,
                                                                       12628, 8266, 12691, 13006,
                                                                       8716, 13069, 5851, 5941,
                                                                       9661, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 13762, 0, 3, 6,
                                                                       12691, 8311, 12754, 13069,
                                                                       8761, 13132, 5941, 6031,
                                                                       9796, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 13951, 0, 6,
                                                                       12817, 8491, 12880, 6211,
                                                                       6271, 10111, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 14077, 0, 6,
                                                                       12880, 8536, 12943, 6271,
                                                                       6331, 10201, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 14203, 0, 6,
                                                                       13006, 8716, 13069, 6451,
                                                                       6511, 10471, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 14329, 0, 6,
                                                                       13069, 8761, 13132, 6511,
                                                                       6571, 10561, ncols, gamma,
                                                                       p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 14455, 0, 3, 6,
                                                                       13195, 9121, 13384, 13951,
                                                                       10111, 14077, 6691, 6871,
                                                                       11191, ncols, gamma, p,
                                                                       q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 14833, 0, 3, 6,
                                                                       13573, 9661, 13762, 14203,
                                                                       10471, 14329, 7231, 7411,
                                                                       12001, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15211, 6, 7771,
                                                                       7786, 12271, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15239, 6, 7786,
                                                                       7801, 12292, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15267, 6, 7801,
                                                                       7816, 12313, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15295, 6, 7816,
                                                                       7831, 12334, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15323, 6, 7861,
                                                                       7876, 12355, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15351, 6, 7876,
                                                                       7891, 12376, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15379, 6, 7891,
                                                                       7906, 12397, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15407, 6, 7906,
                                                                       7921, 12418, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 15435, 3, 6,
                                                                       15211, 12271, 15239, 7951,
                                                                       7996, 12439, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 15519, 3, 6,
                                                                       15239, 12292, 15267, 7996,
                                                                       8041, 12502, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 15603, 3, 6,
                                                                       15267, 12313, 15295, 8041,
                                                                       8086, 12565, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 15687, 3, 6,
                                                                       15323, 12355, 15351, 8176,
                                                                       8221, 12628, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 15771, 3, 6,
                                                                       15351, 12376, 15379, 8221,
                                                                       8266, 12691, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 15855, 3, 6,
                                                                       15379, 12397, 15407, 8266,
                                                                       8311, 12754, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 15939, 0, 6,
                                                                       15211, 12271, 15239, 8401,
                                                                       8446, 12817, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16023, 0, 6,
                                                                       15239, 12292, 15267, 8446,
                                                                       8491, 12880, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16107, 0, 6,
                                                                       15267, 12313, 15295, 8491,
                                                                       8536, 12943, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16191, 0, 6,
                                                                       15323, 12355, 15351, 8626,
                                                                       8671, 13006, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16275, 0, 6,
                                                                       15351, 12376, 15379, 8671,
                                                                       8716, 13069, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16359, 0, 6,
                                                                       15379, 12397, 15407, 8716,
                                                                       8761, 13132, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 16443, 0, 3, 6,
                                                                       15435, 12439, 15519,
                                                                       15939, 12817, 16023, 8851,
                                                                       8986, 13195, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 16695, 0, 3, 6,
                                                                       15519, 12502, 15603,
                                                                       16023, 12880, 16107, 8986,
                                                                       9121, 13384, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 16947, 0, 3, 6,
                                                                       15687, 12628, 15771,
                                                                       16191, 13006, 16275, 9391,
                                                                       9526, 13573, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 17199, 0, 3, 6,
                                                                       15771, 12691, 15855,
                                                                       16275, 13069, 16359, 9526,
                                                                       9661, 13762, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 17451, 0, 6,
                                                                       15939, 12817, 16023, 9931,
                                                                       10021, 13951, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 17619, 0, 6,
                                                                       16023, 12880, 16107,
                                                                       10021, 10111, 14077,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 17787, 0, 6,
                                                                       16191, 13006, 16275,
                                                                       10291, 10381, 14203,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 17955, 0, 6,
                                                                       16275, 13069, 16359,
                                                                       10381, 10471, 14329,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 18123, 0, 3, 6,
                                                                       16443, 13195, 16695,
                                                                       17451, 13951, 17619,
                                                                       10651, 10921, 14455,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 18627, 0, 3, 6,
                                                                       16947, 13573, 17199,
                                                                       17787, 14203, 17955,
                                                                       11461, 11731, 14833,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_s_x(buffer, 19131, 18627, 6, 28, ncols, beta);

                    simdgeo::geom_s_y(buffer, 19299, 18627, 6, 28, ncols, beta);

                    simdgeo::geom_s_z(buffer, 19467, 18627, 6, 28, ncols, beta);

                    simdgeo::geom_s_x(buffer, 19635, 18123, 6, 28, ncols, beta);

                    simdgeo::geom_s_y(buffer, 19803, 18123, 6, 28, ncols, beta);

                    simdgeo::geom_s_z(buffer, 19971, 18123, 6, 28, ncols, beta);

                    simdfunc::contract_primitives(buffer, 20139, 19131, 1008, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 21147, 20139, 6, 1, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 21147, 13, nmax);

        simdtrf::transform_i_inner(buffer, 21147, 20307, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 65 * nvalues + n * npairs, nvalues, buffer, 21147,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 21147, 20475, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 130 * nvalues + n * npairs, nvalues, buffer, 21147,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 21147, 20643, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 195 * nvalues + n * npairs, nvalues, buffer, 21147,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 21147, 20811, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 260 * nvalues + n * npairs, nvalues, buffer, 21147,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 21147, 20979, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 325 * nvalues + n * npairs, nvalues, buffer, 21147,
                                   13, nmax);
    }

    for (size_t m = 0; m < 390; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
