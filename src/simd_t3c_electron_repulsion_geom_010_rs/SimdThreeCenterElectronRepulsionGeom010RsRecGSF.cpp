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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGSF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformG.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_gsf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_gsf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 16720, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 16720, 15715, 900, dimensions);

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
                                                            4, 5, 6, 7, 8}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 18, 6, {1, 2, 3, 4,
                                                        5, 6, 7, 8}, ncols, fj, i * nprim_b + j,
                                                        fq);

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

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1201, 6, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1204, 6, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1207, 6, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1210, 6, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1213, 6, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1216, 6, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1219, 6, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1222, 6, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1225, 6, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1228, 6, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1231, 6, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1234, 6, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1237, 6, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1240, 6, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1243, 6, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1246, 6, 26,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1249, 6, 12, 33,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1258, 6, 13, 36,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1267, 6, 14, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1276, 6, 15, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1285, 6, 16, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1294, 6, 21, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1303, 6, 22, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1312, 6, 23, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1321, 6, 24, 63,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1330, 6, 25, 66,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1339, 6, 10, 69,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1348, 6, 11, 72,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1357, 6, 12, 75,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1366, 6, 13, 78,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1375, 6, 14, 81,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1384, 6, 15, 84,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1393, 6, 16, 87,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1402, 6, 19, 90,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1411, 6, 20, 93,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1420, 6, 21, 96,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1429, 6, 22, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1438, 6, 23, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1447, 6, 24, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1456, 6, 25, 108,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1465, 6, 27, 69,
                                                                       111, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1492, 6, 30, 72,
                                                                       120, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1519, 6, 33, 75,
                                                                       129, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1546, 6, 36, 78,
                                                                       138, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1573, 6, 39, 81,
                                                                       147, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1600, 6, 42, 84,
                                                                       156, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1627, 6, 48, 90,
                                                                       165, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1654, 6, 51, 93,
                                                                       174, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1681, 6, 54, 96,
                                                                       183, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1708, 6, 57, 99,
                                                                       192, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1735, 6, 60, 102,
                                                                       201, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1762, 6, 63, 105,
                                                                       210, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1789, 6, 69, 219,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1807, 6, 72, 225,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1825, 6, 75, 231,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1843, 6, 78, 237,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1861, 6, 81, 243,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1879, 6, 84, 249,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1897, 6, 90, 255,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1915, 6, 93, 261,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1933, 6, 96, 267,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1951, 6, 99, 273,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1969, 6, 102, 279,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1987, 6, 105, 285,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2005, 0, 6, 1465,
                                                                       111, 1492, 219, 291,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2059, 0, 6, 1492,
                                                                       120, 1519, 225, 309,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2113, 0, 6, 1519,
                                                                       129, 1546, 231, 327,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2167, 0, 6, 1546,
                                                                       138, 1573, 237, 345,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2221, 0, 6, 1573,
                                                                       147, 1600, 243, 363,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2275, 0, 6, 1627,
                                                                       165, 1654, 255, 381,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2329, 0, 6, 1654,
                                                                       174, 1681, 261, 399,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2383, 0, 6, 1681,
                                                                       183, 1708, 267, 417,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2437, 0, 6, 1708,
                                                                       192, 1735, 273, 435,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2491, 0, 6, 1735,
                                                                       201, 1762, 279, 453,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2545, 6, 219, 471,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2575, 6, 225, 481,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2605, 6, 231, 491,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2635, 6, 237, 501,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2665, 6, 243, 511,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2695, 6, 255, 521,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2725, 6, 261, 531,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2755, 6, 267, 541,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2785, 6, 273, 551,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2815, 6, 279, 561,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2845, 0, 6, 2005,
                                                                       291, 2059, 471, 571,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2935, 0, 6, 2059,
                                                                       309, 2113, 481, 601,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3025, 0, 6, 2113,
                                                                       327, 2167, 491, 631,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3115, 0, 6, 2167,
                                                                       345, 2221, 501, 661,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3205, 0, 6, 2275,
                                                                       381, 2329, 521, 691,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3295, 0, 6, 2329,
                                                                       399, 2383, 531, 721,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3385, 0, 6, 2383,
                                                                       417, 2437, 541, 751,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3475, 0, 6, 2437,
                                                                       435, 2491, 551, 781,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3565, 6, 471, 811,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3610, 6, 481, 826,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3655, 6, 491, 841,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3700, 6, 501, 856,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3745, 6, 521, 871,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3790, 6, 531, 886,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3835, 6, 541, 901,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3880, 6, 551, 916,
                                                                       ncols, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 3925, 0, 6, 2845,
                                                                       571, 2935, 811, 931,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 4060, 0, 6, 2935,
                                                                       601, 3025, 826, 976,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 4195, 0, 6, 3025,
                                                                       631, 3115, 841, 1021,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 4330, 0, 6, 3205,
                                                                       691, 3295, 871, 1066,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 4465, 0, 6, 3295,
                                                                       721, 3385, 886, 1111,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 4600, 0, 6, 3385,
                                                                       751, 3475, 901, 1156,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4735, 6, 10, 11,
                                                                       1207, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4741, 6, 11, 12,
                                                                       1210, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4747, 6, 12, 13,
                                                                       1213, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4753, 6, 13, 14,
                                                                       1216, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4759, 6, 14, 15,
                                                                       1219, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4765, 6, 15, 16,
                                                                       1222, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4771, 6, 19, 20,
                                                                       1231, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4777, 6, 20, 21,
                                                                       1234, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4783, 6, 21, 22,
                                                                       1237, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4789, 6, 22, 23,
                                                                       1240, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4795, 6, 23, 24,
                                                                       1243, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4801, 6, 24, 25,
                                                                       1246, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4807, 3, 6, 4735,
                                                                       1207, 4741, 1249, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4825, 3, 6, 4741,
                                                                       1210, 4747, 1258, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4843, 3, 6, 4747,
                                                                       1213, 4753, 1267, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4861, 3, 6, 4753,
                                                                       1216, 4759, 1276, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4879, 3, 6, 4759,
                                                                       1219, 4765, 1285, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4897, 3, 6, 4771,
                                                                       1231, 4777, 1294, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4915, 3, 6, 4777,
                                                                       1234, 4783, 1303, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4933, 3, 6, 4783,
                                                                       1237, 4789, 1312, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4951, 3, 6, 4789,
                                                                       1240, 4795, 1321, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4969, 3, 6, 4795,
                                                                       1243, 4801, 1330, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4987, 0, 6, 4735,
                                                                       1207, 4741, 1357, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5005, 0, 6, 4741,
                                                                       1210, 4747, 1366, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5023, 0, 6, 4747,
                                                                       1213, 4753, 1375, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5041, 0, 6, 4753,
                                                                       1216, 4759, 1384, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5059, 0, 6, 4759,
                                                                       1219, 4765, 1393, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5077, 0, 6, 4771,
                                                                       1231, 4777, 1420, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5095, 0, 6, 4777,
                                                                       1234, 4783, 1429, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5113, 0, 6, 4783,
                                                                       1237, 4789, 1438, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5131, 0, 6, 4789,
                                                                       1240, 4795, 1447, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5149, 0, 6, 4795,
                                                                       1243, 4801, 1456, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5167, 0, 3, 6,
                                                                       4807, 1249, 4825, 4987,
                                                                       1357, 5005, 111, 120,
                                                                       1519, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5221, 0, 3, 6,
                                                                       4825, 1258, 4843, 5005,
                                                                       1366, 5023, 120, 129,
                                                                       1546, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5275, 0, 3, 6,
                                                                       4843, 1267, 4861, 5023,
                                                                       1375, 5041, 129, 138,
                                                                       1573, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5329, 0, 3, 6,
                                                                       4861, 1276, 4879, 5041,
                                                                       1384, 5059, 138, 147,
                                                                       1600, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5383, 0, 3, 6,
                                                                       4897, 1294, 4915, 5077,
                                                                       1420, 5095, 165, 174,
                                                                       1681, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5437, 0, 3, 6,
                                                                       4915, 1303, 4933, 5095,
                                                                       1429, 5113, 174, 183,
                                                                       1708, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5491, 0, 3, 6,
                                                                       4933, 1312, 4951, 5113,
                                                                       1438, 5131, 183, 192,
                                                                       1735, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5545, 0, 3, 6,
                                                                       4951, 1321, 4969, 5131,
                                                                       1447, 5149, 192, 201,
                                                                       1762, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5599, 0, 6, 4987,
                                                                       1357, 5005, 219, 225,
                                                                       1825, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5635, 0, 6, 5005,
                                                                       1366, 5023, 225, 231,
                                                                       1843, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5671, 0, 6, 5023,
                                                                       1375, 5041, 231, 237,
                                                                       1861, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5707, 0, 6, 5041,
                                                                       1384, 5059, 237, 243,
                                                                       1879, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5743, 0, 6, 5077,
                                                                       1420, 5095, 255, 261,
                                                                       1933, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5779, 0, 6, 5095,
                                                                       1429, 5113, 261, 267,
                                                                       1951, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5815, 0, 6, 5113,
                                                                       1438, 5131, 267, 273,
                                                                       1969, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5851, 0, 6, 5131,
                                                                       1447, 5149, 273, 279,
                                                                       1987, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5887, 0, 3, 6,
                                                                       5167, 1519, 5221, 5599,
                                                                       1825, 5635, 291, 309,
                                                                       2113, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5995, 0, 3, 6,
                                                                       5221, 1546, 5275, 5635,
                                                                       1843, 5671, 309, 327,
                                                                       2167, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 6103, 0, 3, 6,
                                                                       5275, 1573, 5329, 5671,
                                                                       1861, 5707, 327, 345,
                                                                       2221, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 6211, 0, 3, 6,
                                                                       5383, 1681, 5437, 5743,
                                                                       1933, 5779, 381, 399,
                                                                       2383, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 6319, 0, 3, 6,
                                                                       5437, 1708, 5491, 5779,
                                                                       1951, 5815, 399, 417,
                                                                       2437, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 6427, 0, 3, 6,
                                                                       5491, 1735, 5545, 5815,
                                                                       1969, 5851, 417, 435,
                                                                       2491, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6535, 0, 6, 5599,
                                                                       1825, 5635, 471, 481,
                                                                       2605, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6595, 0, 6, 5635,
                                                                       1843, 5671, 481, 491,
                                                                       2635, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6655, 0, 6, 5671,
                                                                       1861, 5707, 491, 501,
                                                                       2665, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6715, 0, 6, 5743,
                                                                       1933, 5779, 521, 531,
                                                                       2755, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6775, 0, 6, 5779,
                                                                       1951, 5815, 531, 541,
                                                                       2785, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6835, 0, 6, 5815,
                                                                       1969, 5851, 541, 551,
                                                                       2815, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6895, 0, 3, 6,
                                                                       5167, 5221, 5887, 2113,
                                                                       5995, 6535, 2605, 6595,
                                                                       571, 601, 3025, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 7075, 0, 3, 6,
                                                                       5221, 5275, 5995, 2167,
                                                                       6103, 6595, 2635, 6655,
                                                                       601, 631, 3115, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 7255, 0, 3, 6,
                                                                       5383, 5437, 6211, 2383,
                                                                       6319, 6715, 2755, 6775,
                                                                       691, 721, 3385, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 7435, 0, 3, 6,
                                                                       5437, 5491, 6319, 2437,
                                                                       6427, 6775, 2785, 6835,
                                                                       721, 751, 3475, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7615, 0, 6, 6535,
                                                                       2605, 6595, 811, 826,
                                                                       3655, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7705, 0, 6, 6595,
                                                                       2635, 6655, 826, 841,
                                                                       3700, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7795, 0, 6, 6715,
                                                                       2755, 6775, 871, 886,
                                                                       3835, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7885, 0, 6, 6775,
                                                                       2785, 6835, 886, 901,
                                                                       3880, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 7975, 0, 3, 6,
                                                                       5887, 5995, 6895, 3025,
                                                                       7075, 7615, 3655, 7705,
                                                                       931, 976, 4195, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 8245, 0, 3, 6,
                                                                       6211, 6319, 7255, 3385,
                                                                       7435, 7795, 3835, 7885,
                                                                       1066, 1111, 4600, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8515, 6, 1201,
                                                                       1204, 4735, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8525, 6, 1204,
                                                                       1207, 4741, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8535, 6, 1207,
                                                                       1210, 4747, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8545, 6, 1210,
                                                                       1213, 4753, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8555, 6, 1213,
                                                                       1216, 4759, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8565, 6, 1216,
                                                                       1219, 4765, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8575, 6, 1225,
                                                                       1228, 4771, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8585, 6, 1228,
                                                                       1231, 4777, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8595, 6, 1231,
                                                                       1234, 4783, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8605, 6, 1234,
                                                                       1237, 4789, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8615, 6, 1237,
                                                                       1240, 4795, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 8625, 6, 1240,
                                                                       1243, 4801, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8635, 3, 6, 8515,
                                                                       4735, 8525, 4807, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8665, 3, 6, 8525,
                                                                       4741, 8535, 4825, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8695, 3, 6, 8535,
                                                                       4747, 8545, 4843, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8725, 3, 6, 8545,
                                                                       4753, 8555, 4861, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8755, 3, 6, 8555,
                                                                       4759, 8565, 4879, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8785, 3, 6, 8575,
                                                                       4771, 8585, 4897, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8815, 3, 6, 8585,
                                                                       4777, 8595, 4915, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8845, 3, 6, 8595,
                                                                       4783, 8605, 4933, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8875, 3, 6, 8605,
                                                                       4789, 8615, 4951, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 8905, 3, 6, 8615,
                                                                       4795, 8625, 4969, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8935, 0, 6, 8515,
                                                                       4735, 8525, 1339, 1348,
                                                                       4987, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8965, 0, 6, 8525,
                                                                       4741, 8535, 1348, 1357,
                                                                       5005, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8995, 0, 6, 8535,
                                                                       4747, 8545, 1357, 1366,
                                                                       5023, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9025, 0, 6, 8545,
                                                                       4753, 8555, 1366, 1375,
                                                                       5041, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9055, 0, 6, 8555,
                                                                       4759, 8565, 1375, 1384,
                                                                       5059, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9085, 0, 6, 8575,
                                                                       4771, 8585, 1402, 1411,
                                                                       5077, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9115, 0, 6, 8585,
                                                                       4777, 8595, 1411, 1420,
                                                                       5095, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9145, 0, 6, 8595,
                                                                       4783, 8605, 1420, 1429,
                                                                       5113, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9175, 0, 6, 8605,
                                                                       4789, 8615, 1429, 1438,
                                                                       5131, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9205, 0, 6, 8615,
                                                                       4795, 8625, 1438, 1447,
                                                                       5149, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9235, 0, 3, 6,
                                                                       8635, 4807, 8665, 8935,
                                                                       4987, 8965, 1465, 1492,
                                                                       5167, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9325, 0, 3, 6,
                                                                       8665, 4825, 8695, 8965,
                                                                       5005, 8995, 1492, 1519,
                                                                       5221, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9415, 0, 3, 6,
                                                                       8695, 4843, 8725, 8995,
                                                                       5023, 9025, 1519, 1546,
                                                                       5275, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9505, 0, 3, 6,
                                                                       8725, 4861, 8755, 9025,
                                                                       5041, 9055, 1546, 1573,
                                                                       5329, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9595, 0, 3, 6,
                                                                       8785, 4897, 8815, 9085,
                                                                       5077, 9115, 1627, 1654,
                                                                       5383, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9685, 0, 3, 6,
                                                                       8815, 4915, 8845, 9115,
                                                                       5095, 9145, 1654, 1681,
                                                                       5437, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9775, 0, 3, 6,
                                                                       8845, 4933, 8875, 9145,
                                                                       5113, 9175, 1681, 1708,
                                                                       5491, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9865, 0, 3, 6,
                                                                       8875, 4951, 8905, 9175,
                                                                       5131, 9205, 1708, 1735,
                                                                       5545, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 9955, 0, 6, 8935,
                                                                       4987, 8965, 1789, 1807,
                                                                       5599, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10015, 0, 6, 8965,
                                                                       5005, 8995, 1807, 1825,
                                                                       5635, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10075, 0, 6, 8995,
                                                                       5023, 9025, 1825, 1843,
                                                                       5671, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10135, 0, 6, 9025,
                                                                       5041, 9055, 1843, 1861,
                                                                       5707, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10195, 0, 6, 9085,
                                                                       5077, 9115, 1897, 1915,
                                                                       5743, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10255, 0, 6, 9115,
                                                                       5095, 9145, 1915, 1933,
                                                                       5779, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10315, 0, 6, 9145,
                                                                       5113, 9175, 1933, 1951,
                                                                       5815, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10375, 0, 6, 9175,
                                                                       5131, 9205, 1951, 1969,
                                                                       5851, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10435, 0, 3, 6,
                                                                       9235, 5167, 9325, 9955,
                                                                       5599, 10015, 2005, 2059,
                                                                       5887, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10615, 0, 3, 6,
                                                                       9325, 5221, 9415, 10015,
                                                                       5635, 10075, 2059, 2113,
                                                                       5995, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10795, 0, 3, 6,
                                                                       9415, 5275, 9505, 10075,
                                                                       5671, 10135, 2113, 2167,
                                                                       6103, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10975, 0, 3, 6,
                                                                       9595, 5383, 9685, 10195,
                                                                       5743, 10255, 2275, 2329,
                                                                       6211, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 11155, 0, 3, 6,
                                                                       9685, 5437, 9775, 10255,
                                                                       5779, 10315, 2329, 2383,
                                                                       6319, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 11335, 0, 3, 6,
                                                                       9775, 5491, 9865, 10315,
                                                                       5815, 10375, 2383, 2437,
                                                                       6427, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11515, 0, 6, 9955,
                                                                       5599, 10015, 2545, 2575,
                                                                       6535, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11615, 0, 6,
                                                                       10015, 5635, 10075, 2575,
                                                                       2605, 6595, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11715, 0, 6,
                                                                       10075, 5671, 10135, 2605,
                                                                       2635, 6655, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11815, 0, 6,
                                                                       10195, 5743, 10255, 2695,
                                                                       2725, 6715, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11915, 0, 6,
                                                                       10255, 5779, 10315, 2725,
                                                                       2755, 6775, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 12015, 0, 6,
                                                                       10315, 5815, 10375, 2755,
                                                                       2785, 6835, ncols, gamma,
                                                                       p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 12115, 0, 3, 6,
                                                                       9235, 9325, 10435, 5887,
                                                                       10615, 11515, 6535, 11615,
                                                                       2845, 2935, 6895, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 12415, 0, 3, 6,
                                                                       9325, 9415, 10615, 5995,
                                                                       10795, 11615, 6595, 11715,
                                                                       2935, 3025, 7075, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 12715, 0, 3, 6,
                                                                       9595, 9685, 10975, 6211,
                                                                       11155, 11815, 6715, 11915,
                                                                       3205, 3295, 7255, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 13015, 0, 3, 6,
                                                                       9685, 9775, 11155, 6319,
                                                                       11335, 11915, 6775, 12015,
                                                                       3295, 3385, 7435, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13315, 0, 6,
                                                                       11515, 6535, 11615, 3565,
                                                                       3610, 7615, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13465, 0, 6,
                                                                       11615, 6595, 11715, 3610,
                                                                       3655, 7705, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13615, 0, 6,
                                                                       11815, 6715, 11915, 3745,
                                                                       3790, 7795, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13765, 0, 6,
                                                                       11915, 6775, 12015, 3790,
                                                                       3835, 7885, ncols, gamma,
                                                                       p, q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 13915, 0, 3, 6,
                                                                       10435, 10615, 12115, 6895,
                                                                       12415, 13315, 7615, 13465,
                                                                       3925, 4060, 7975, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 14365, 0, 3, 6,
                                                                       10975, 11155, 12715, 7255,
                                                                       13015, 13615, 7795, 13765,
                                                                       4330, 4465, 8245, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 14815, 14365, 15, 10, ncols, beta);

                    simdgeo::geom_s_y(buffer, 14965, 14365, 15, 10, ncols, beta);

                    simdgeo::geom_s_z(buffer, 15115, 14365, 15, 10, ncols, beta);

                    simdgeo::geom_s_x(buffer, 15265, 13915, 15, 10, ncols, beta);

                    simdgeo::geom_s_y(buffer, 15415, 13915, 15, 10, ncols, beta);

                    simdgeo::geom_s_z(buffer, 15565, 13915, 15, 10, ncols, beta);

                    simdfunc::contract_primitives(buffer, 15715, 14815, 900, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 16615, 15715, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 16615, 7, nmax);

        simdtrf::transform_f_inner(buffer, 16615, 15865, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 16615, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 16615, 16015, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 126 * nvalues + n * npairs, nvalues, buffer, 16615,
                                   7, nmax);

        simdtrf::transform_f_inner(buffer, 16615, 16165, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 189 * nvalues + n * npairs, nvalues, buffer, 16615,
                                   7, nmax);

        simdtrf::transform_f_inner(buffer, 16615, 16315, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 252 * nvalues + n * npairs, nvalues, buffer, 16615,
                                   7, nmax);

        simdtrf::transform_f_inner(buffer, 16615, 16465, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 16615,
                                   7, nmax);
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
