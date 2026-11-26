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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGSD.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformD.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_gsd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_gsd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 8026, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 270 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 8026, 7411, 540, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 9, 6, 7,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 18, 6, 7,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 27, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 30, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 3, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 3, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 3, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 3, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 60, 3, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 63, 3, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 66, 3, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 69, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 72, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 75, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 78, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 81, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 84, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 87, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 90, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 93, 0, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 96, 0, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 99, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 102, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 105, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 108, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 111, 0, 6, 10, 11,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 120, 0, 6, 11, 12,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 129, 0, 6, 12, 13,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 138, 0, 6, 13, 14,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 147, 0, 6, 14, 15,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 156, 0, 6, 15, 16,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 165, 0, 6, 19, 20,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 174, 0, 6, 20, 21,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 183, 0, 6, 21, 22,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 192, 0, 6, 22, 23,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 201, 0, 6, 23, 24,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 210, 0, 6, 24, 25,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 219, 0, 6, 10, 11,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 225, 0, 6, 11, 12,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 231, 0, 6, 12, 13,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 237, 0, 6, 13, 14,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 243, 0, 6, 14, 15,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 249, 0, 6, 15, 16,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 255, 0, 6, 19, 20,
                                                                       90, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 261, 0, 6, 20, 21,
                                                                       93, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 267, 0, 6, 21, 22,
                                                                       96, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 273, 0, 6, 22, 23,
                                                                       99, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 279, 0, 6, 23, 24,
                                                                       102, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 285, 0, 6, 24, 25,
                                                                       105, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 291, 0, 3, 6, 69,
                                                                       72, 111, 120, 219, 225,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 309, 0, 3, 6, 72,
                                                                       75, 120, 129, 225, 231,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 327, 0, 3, 6, 75,
                                                                       78, 129, 138, 231, 237,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 345, 0, 3, 6, 78,
                                                                       81, 138, 147, 237, 243,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 363, 0, 3, 6, 81,
                                                                       84, 147, 156, 243, 249,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 381, 0, 3, 6, 90,
                                                                       93, 165, 174, 255, 261,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 399, 0, 3, 6, 93,
                                                                       96, 174, 183, 261, 267,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 417, 0, 3, 6, 96,
                                                                       99, 183, 192, 267, 273,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 435, 0, 3, 6, 99,
                                                                       102, 192, 201, 273, 279,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 453, 0, 3, 6, 102,
                                                                       105, 201, 210, 279, 285,
                                                                       ncols, gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 471, 0, 6, 69, 72,
                                                                       219, 225, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 481, 0, 6, 72, 75,
                                                                       225, 231, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 491, 0, 6, 75, 78,
                                                                       231, 237, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 501, 0, 6, 78, 81,
                                                                       237, 243, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 511, 0, 6, 81, 84,
                                                                       243, 249, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 521, 0, 6, 90, 93,
                                                                       255, 261, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 531, 0, 6, 93, 96,
                                                                       261, 267, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 541, 0, 6, 96, 99,
                                                                       267, 273, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 551, 0, 6, 99,
                                                                       102, 273, 279, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 561, 0, 6, 102,
                                                                       105, 279, 285, ncols,
                                                                       gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 571, 0, 3, 6, 219,
                                                                       225, 291, 309, 471, 481,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 601, 0, 3, 6, 225,
                                                                       231, 309, 327, 481, 491,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 631, 0, 3, 6, 231,
                                                                       237, 327, 345, 491, 501,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 661, 0, 3, 6, 237,
                                                                       243, 345, 363, 501, 511,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 691, 0, 3, 6, 255,
                                                                       261, 381, 399, 521, 531,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 721, 0, 3, 6, 261,
                                                                       267, 399, 417, 531, 541,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 751, 0, 3, 6, 267,
                                                                       273, 417, 435, 541, 551,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 781, 0, 3, 6, 273,
                                                                       279, 435, 453, 551, 561,
                                                                       ncols, gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 811, 0, 6, 219,
                                                                       225, 471, 481, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 826, 0, 6, 225,
                                                                       231, 481, 491, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 841, 0, 6, 231,
                                                                       237, 491, 501, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 856, 0, 6, 237,
                                                                       243, 501, 511, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 871, 0, 6, 255,
                                                                       261, 521, 531, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 886, 0, 6, 261,
                                                                       267, 531, 541, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 901, 0, 6, 267,
                                                                       273, 541, 551, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 916, 0, 6, 273,
                                                                       279, 551, 561, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 931, 0, 3, 6, 291,
                                                                       309, 471, 481, 571, 601,
                                                                       811, 826, ncols, gamma, p,
                                                                       q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 976, 0, 3, 6, 309,
                                                                       327, 481, 491, 601, 631,
                                                                       826, 841, ncols, gamma, p,
                                                                       q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1021, 0, 3, 6,
                                                                       327, 345, 491, 501, 631,
                                                                       661, 841, 856, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1066, 0, 3, 6,
                                                                       381, 399, 521, 531, 691,
                                                                       721, 871, 886, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1111, 0, 3, 6,
                                                                       399, 417, 531, 541, 721,
                                                                       751, 886, 901, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 6,
                                                                       417, 435, 541, 551, 751,
                                                                       781, 901, 916, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1201, 6, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1204, 6, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1207, 6, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1210, 6, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1213, 6, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1216, 6, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1219, 6, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1222, 6, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1225, 6, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1228, 6, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1231, 6, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1234, 6, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1237, 6, 12, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1246, 6, 13, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1255, 6, 14, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1264, 6, 15, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1273, 6, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1282, 6, 21, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1291, 6, 22, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1300, 6, 23, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1309, 6, 24, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1318, 6, 25, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1327, 6, 12, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1336, 6, 13, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1345, 6, 14, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1354, 6, 15, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1363, 6, 16, 87,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1372, 6, 21, 96,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1381, 6, 22, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1390, 6, 23, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1399, 6, 24, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1408, 6, 25, 108,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1417, 6, 33, 75,
                                                                       129, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1444, 6, 36, 78,
                                                                       138, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1471, 6, 39, 81,
                                                                       147, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1498, 6, 42, 84,
                                                                       156, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1525, 6, 54, 96,
                                                                       183, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1552, 6, 57, 99,
                                                                       192, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1579, 6, 60, 102,
                                                                       201, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1606, 6, 63, 105,
                                                                       210, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1633, 6, 75, 231,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1651, 6, 78, 237,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1669, 6, 81, 243,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1687, 6, 84, 249,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1705, 6, 96, 267,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1723, 6, 99, 273,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1741, 6, 102, 279,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1759, 6, 105, 285,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1777, 0, 6, 1417,
                                                                       129, 1444, 231, 327,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1831, 0, 6, 1444,
                                                                       138, 1471, 237, 345,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1885, 0, 6, 1471,
                                                                       147, 1498, 243, 363,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1939, 0, 6, 1525,
                                                                       183, 1552, 267, 417,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1993, 0, 6, 1552,
                                                                       192, 1579, 273, 435,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2047, 0, 6, 1579,
                                                                       201, 1606, 279, 453,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2101, 6, 231, 491,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2131, 6, 237, 501,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2161, 6, 243, 511,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2191, 6, 267, 541,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2221, 6, 273, 551,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2251, 6, 279, 561,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2281, 0, 6, 1777,
                                                                       327, 1831, 491, 631,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2371, 0, 6, 1831,
                                                                       345, 1885, 501, 661,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2461, 0, 6, 1939,
                                                                       417, 1993, 541, 751,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2551, 0, 6, 1993,
                                                                       435, 2047, 551, 781,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2641, 6, 491, 841,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2686, 6, 501, 856,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2731, 6, 541, 901,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2776, 6, 551, 916,
                                                                       ncols, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 2821, 0, 6, 2281,
                                                                       631, 2371, 841, 1021,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 2956, 0, 6, 2461,
                                                                       751, 2551, 901, 1156,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3091, 6, 10, 11,
                                                                       1201, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3097, 6, 11, 12,
                                                                       1204, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3103, 6, 12, 13,
                                                                       1207, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3109, 6, 13, 14,
                                                                       1210, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3115, 6, 14, 15,
                                                                       1213, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3121, 6, 15, 16,
                                                                       1216, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3127, 6, 19, 20,
                                                                       1219, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3133, 6, 20, 21,
                                                                       1222, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3139, 6, 21, 22,
                                                                       1225, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3145, 6, 22, 23,
                                                                       1228, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3151, 6, 23, 24,
                                                                       1231, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3157, 6, 24, 25,
                                                                       1234, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3163, 3, 6, 3091,
                                                                       1201, 3097, 1237, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3181, 3, 6, 3097,
                                                                       1204, 3103, 1246, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3199, 3, 6, 3103,
                                                                       1207, 3109, 1255, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3217, 3, 6, 3109,
                                                                       1210, 3115, 1264, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3235, 3, 6, 3115,
                                                                       1213, 3121, 1273, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3253, 3, 6, 3127,
                                                                       1219, 3133, 1282, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3271, 3, 6, 3133,
                                                                       1222, 3139, 1291, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3289, 3, 6, 3139,
                                                                       1225, 3145, 1300, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3307, 3, 6, 3145,
                                                                       1228, 3151, 1309, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3325, 3, 6, 3151,
                                                                       1231, 3157, 1318, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3343, 0, 6, 3091,
                                                                       1201, 3097, 1327, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3361, 0, 6, 3097,
                                                                       1204, 3103, 1336, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3379, 0, 6, 3103,
                                                                       1207, 3109, 1345, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3397, 0, 6, 3109,
                                                                       1210, 3115, 1354, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3415, 0, 6, 3115,
                                                                       1213, 3121, 1363, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3433, 0, 6, 3127,
                                                                       1219, 3133, 1372, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3451, 0, 6, 3133,
                                                                       1222, 3139, 1381, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3469, 0, 6, 3139,
                                                                       1225, 3145, 1390, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3487, 0, 6, 3145,
                                                                       1228, 3151, 1399, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3505, 0, 6, 3151,
                                                                       1231, 3157, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3523, 0, 3, 6,
                                                                       3163, 1237, 3181, 3343,
                                                                       1327, 3361, 111, 120,
                                                                       1417, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3577, 0, 3, 6,
                                                                       3181, 1246, 3199, 3361,
                                                                       1336, 3379, 120, 129,
                                                                       1444, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3631, 0, 3, 6,
                                                                       3199, 1255, 3217, 3379,
                                                                       1345, 3397, 129, 138,
                                                                       1471, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3685, 0, 3, 6,
                                                                       3217, 1264, 3235, 3397,
                                                                       1354, 3415, 138, 147,
                                                                       1498, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3739, 0, 3, 6,
                                                                       3253, 1282, 3271, 3433,
                                                                       1372, 3451, 165, 174,
                                                                       1525, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3793, 0, 3, 6,
                                                                       3271, 1291, 3289, 3451,
                                                                       1381, 3469, 174, 183,
                                                                       1552, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3847, 0, 3, 6,
                                                                       3289, 1300, 3307, 3469,
                                                                       1390, 3487, 183, 192,
                                                                       1579, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3901, 0, 3, 6,
                                                                       3307, 1309, 3325, 3487,
                                                                       1399, 3505, 192, 201,
                                                                       1606, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3955, 0, 6, 3343,
                                                                       1327, 3361, 219, 225,
                                                                       1633, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3991, 0, 6, 3361,
                                                                       1336, 3379, 225, 231,
                                                                       1651, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4027, 0, 6, 3379,
                                                                       1345, 3397, 231, 237,
                                                                       1669, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4063, 0, 6, 3397,
                                                                       1354, 3415, 237, 243,
                                                                       1687, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4099, 0, 6, 3433,
                                                                       1372, 3451, 255, 261,
                                                                       1705, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4135, 0, 6, 3451,
                                                                       1381, 3469, 261, 267,
                                                                       1723, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4171, 0, 6, 3469,
                                                                       1390, 3487, 267, 273,
                                                                       1741, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4207, 0, 6, 3487,
                                                                       1399, 3505, 273, 279,
                                                                       1759, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4243, 0, 3, 6,
                                                                       3523, 1417, 3577, 3955,
                                                                       1633, 3991, 291, 309,
                                                                       1777, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4351, 0, 3, 6,
                                                                       3577, 1444, 3631, 3991,
                                                                       1651, 4027, 309, 327,
                                                                       1831, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4459, 0, 3, 6,
                                                                       3631, 1471, 3685, 4027,
                                                                       1669, 4063, 327, 345,
                                                                       1885, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4567, 0, 3, 6,
                                                                       3739, 1525, 3793, 4099,
                                                                       1705, 4135, 381, 399,
                                                                       1939, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4675, 0, 3, 6,
                                                                       3793, 1552, 3847, 4135,
                                                                       1723, 4171, 399, 417,
                                                                       1993, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4783, 0, 3, 6,
                                                                       3847, 1579, 3901, 4171,
                                                                       1741, 4207, 417, 435,
                                                                       2047, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4891, 0, 6, 3955,
                                                                       1633, 3991, 471, 481,
                                                                       2101, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4951, 0, 6, 3991,
                                                                       1651, 4027, 481, 491,
                                                                       2131, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5011, 0, 6, 4027,
                                                                       1669, 4063, 491, 501,
                                                                       2161, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5071, 0, 6, 4099,
                                                                       1705, 4135, 521, 531,
                                                                       2191, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5131, 0, 6, 4135,
                                                                       1723, 4171, 531, 541,
                                                                       2221, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5191, 0, 6, 4171,
                                                                       1741, 4207, 541, 551,
                                                                       2251, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 5251, 0, 3, 6,
                                                                       3523, 3577, 4243, 1777,
                                                                       4351, 4891, 2101, 4951,
                                                                       571, 601, 2281, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 5431, 0, 3, 6,
                                                                       3577, 3631, 4351, 1831,
                                                                       4459, 4951, 2131, 5011,
                                                                       601, 631, 2371, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 5611, 0, 3, 6,
                                                                       3739, 3793, 4567, 1939,
                                                                       4675, 5071, 2191, 5131,
                                                                       691, 721, 2461, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 5791, 0, 3, 6,
                                                                       3793, 3847, 4675, 1993,
                                                                       4783, 5131, 2221, 5191,
                                                                       721, 751, 2551, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 5971, 0, 6, 4891,
                                                                       2101, 4951, 811, 826,
                                                                       2641, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6061, 0, 6, 4951,
                                                                       2131, 5011, 826, 841,
                                                                       2686, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6151, 0, 6, 5071,
                                                                       2191, 5131, 871, 886,
                                                                       2731, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 6241, 0, 6, 5131,
                                                                       2221, 5191, 886, 901,
                                                                       2776, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 6331, 0, 3, 6,
                                                                       4243, 4351, 5251, 2281,
                                                                       5431, 5971, 2641, 6061,
                                                                       931, 976, 2821, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 6601, 0, 3, 6,
                                                                       4567, 4675, 5611, 2461,
                                                                       5791, 6151, 2731, 6241,
                                                                       1066, 1111, 2956, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 6871, 6601, 15, 6, ncols, beta);

                    simdgeo::geom_s_y(buffer, 6961, 6601, 15, 6, ncols, beta);

                    simdgeo::geom_s_z(buffer, 7051, 6601, 15, 6, ncols, beta);

                    simdgeo::geom_s_x(buffer, 7141, 6331, 15, 6, ncols, beta);

                    simdgeo::geom_s_y(buffer, 7231, 6331, 15, 6, ncols, beta);

                    simdgeo::geom_s_z(buffer, 7321, 6331, 15, 6, ncols, beta);

                    simdfunc::contract_primitives(buffer, 7411, 6871, 540, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 7951, 7411, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 7951, 5, nmax);

        simdtrf::transform_d_inner(buffer, 7951, 7501, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 45 * nvalues + n * npairs, nvalues, buffer, 7951, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 7951, 7591, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 90 * nvalues + n * npairs, nvalues, buffer, 7951, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 7951, 7681, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 135 * nvalues + n * npairs, nvalues, buffer, 7951, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 7951, 7771, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 180 * nvalues + n * npairs, nvalues, buffer, 7951, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 7951, 7861, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 225 * nvalues + n * npairs, nvalues, buffer, 7951, 5,
                                   nmax);
    }

    for (size_t m = 0; m < 270; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
