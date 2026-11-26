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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecFSF.hpp"

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

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_fsf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_fsf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 8925, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 294 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 8925, 8255, 600, dimensions);

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
                                                            4, 5, 6, 7}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 17, 6, {1, 2, 3, 4,
                                                        5, 6, 7}, ncols, fj, i * nprim_b + j,
                                                        fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 3, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 3, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 3, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 55, 3, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 58, 3, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 67, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 70, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 73, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 76, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 79, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 82, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 85, 0, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 88, 0, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 91, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 94, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 97, 0, 6, 10, 11,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 106, 0, 6, 11, 12,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 115, 0, 6, 12, 13,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 124, 0, 6, 13, 14,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 133, 0, 6, 14, 15,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 142, 0, 6, 18, 19,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 151, 0, 6, 19, 20,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 160, 0, 6, 20, 21,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 169, 0, 6, 21, 22,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 178, 0, 6, 22, 23,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 187, 0, 6, 10, 11,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 193, 0, 6, 11, 12,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 199, 0, 6, 12, 13,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 205, 0, 6, 13, 14,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 211, 0, 6, 14, 15,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 217, 0, 6, 18, 19,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 223, 0, 6, 19, 20,
                                                                       82, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 229, 0, 6, 20, 21,
                                                                       85, 88, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 235, 0, 6, 21, 22,
                                                                       88, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 241, 0, 6, 22, 23,
                                                                       91, 94, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 247, 0, 3, 6, 61,
                                                                       64, 97, 106, 187, 193,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 265, 0, 3, 6, 64,
                                                                       67, 106, 115, 193, 199,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 283, 0, 3, 6, 67,
                                                                       70, 115, 124, 199, 205,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 301, 0, 3, 6, 70,
                                                                       73, 124, 133, 205, 211,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 319, 0, 3, 6, 79,
                                                                       82, 142, 151, 217, 223,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 337, 0, 3, 6, 82,
                                                                       85, 151, 160, 223, 229,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 355, 0, 3, 6, 85,
                                                                       88, 160, 169, 229, 235,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 373, 0, 3, 6, 88,
                                                                       91, 169, 178, 235, 241,
                                                                       ncols, gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 391, 0, 6, 61, 64,
                                                                       187, 193, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 401, 0, 6, 64, 67,
                                                                       193, 199, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 411, 0, 6, 67, 70,
                                                                       199, 205, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 421, 0, 6, 70, 73,
                                                                       205, 211, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 431, 0, 6, 79, 82,
                                                                       217, 223, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 441, 0, 6, 82, 85,
                                                                       223, 229, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 451, 0, 6, 85, 88,
                                                                       229, 235, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 461, 0, 6, 88, 91,
                                                                       235, 241, ncols, gamma, p,
                                                                       q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 471, 0, 3, 6, 187,
                                                                       193, 247, 265, 391, 401,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 501, 0, 3, 6, 193,
                                                                       199, 265, 283, 401, 411,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 531, 0, 3, 6, 199,
                                                                       205, 283, 301, 411, 421,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 561, 0, 3, 6, 217,
                                                                       223, 319, 337, 431, 441,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 591, 0, 3, 6, 223,
                                                                       229, 337, 355, 441, 451,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 621, 0, 3, 6, 229,
                                                                       235, 355, 373, 451, 461,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 651, 6, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 654, 6, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 657, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 660, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 663, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 666, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 669, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 672, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 675, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 678, 6, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 681, 6, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 684, 6, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 687, 6, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 690, 6, 24, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 693, 6, 12, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 702, 6, 13, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 711, 6, 14, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 720, 6, 15, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 729, 6, 20, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 738, 6, 21, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 747, 6, 22, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 756, 6, 23, 58,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 765, 6, 10, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 774, 6, 11, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 783, 6, 12, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 792, 6, 13, 70,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 801, 6, 14, 73,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 810, 6, 15, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 819, 6, 18, 79,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 828, 6, 19, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 837, 6, 20, 85,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 846, 6, 21, 88,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 855, 6, 22, 91,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 864, 6, 23, 94,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 873, 6, 25, 61,
                                                                       97, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 900, 6, 28, 64,
                                                                       106, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 927, 6, 31, 67,
                                                                       115, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 954, 6, 34, 70,
                                                                       124, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 981, 6, 37, 73,
                                                                       133, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1008, 6, 43, 79,
                                                                       142, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1035, 6, 46, 82,
                                                                       151, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1062, 6, 49, 85,
                                                                       160, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1089, 6, 52, 88,
                                                                       169, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1116, 6, 55, 91,
                                                                       178, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1143, 6, 61, 187,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1161, 6, 64, 193,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1179, 6, 67, 199,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1197, 6, 70, 205,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1215, 6, 73, 211,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1233, 6, 79, 217,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1251, 6, 82, 223,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1269, 6, 85, 229,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1287, 6, 88, 235,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1305, 6, 91, 241,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1323, 0, 6, 873,
                                                                       97, 900, 187, 247, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1377, 0, 6, 900,
                                                                       106, 927, 193, 265, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1431, 0, 6, 927,
                                                                       115, 954, 199, 283, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1485, 0, 6, 954,
                                                                       124, 981, 205, 301, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1539, 0, 6, 1008,
                                                                       142, 1035, 217, 319,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1593, 0, 6, 1035,
                                                                       151, 1062, 223, 337,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1647, 0, 6, 1062,
                                                                       160, 1089, 229, 355,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1701, 0, 6, 1089,
                                                                       169, 1116, 235, 373,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1755, 6, 187, 391,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1785, 6, 193, 401,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1815, 6, 199, 411,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1845, 6, 205, 421,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1875, 6, 217, 431,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1905, 6, 223, 441,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1935, 6, 229, 451,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 1965, 6, 235, 461,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 1995, 0, 6, 1323,
                                                                       247, 1377, 391, 471,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2085, 0, 6, 1377,
                                                                       265, 1431, 401, 501,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2175, 0, 6, 1431,
                                                                       283, 1485, 411, 531,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2265, 0, 6, 1539,
                                                                       319, 1593, 431, 561,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2355, 0, 6, 1593,
                                                                       337, 1647, 441, 591,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2445, 0, 6, 1647,
                                                                       355, 1701, 451, 621,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2535, 6, 10, 11,
                                                                       657, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2541, 6, 11, 12,
                                                                       660, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2547, 6, 12, 13,
                                                                       663, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2553, 6, 13, 14,
                                                                       666, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2559, 6, 14, 15,
                                                                       669, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2565, 6, 18, 19,
                                                                       678, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2571, 6, 19, 20,
                                                                       681, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2577, 6, 20, 21,
                                                                       684, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2583, 6, 21, 22,
                                                                       687, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2589, 6, 22, 23,
                                                                       690, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2595, 3, 6, 2535,
                                                                       657, 2541, 693, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2613, 3, 6, 2541,
                                                                       660, 2547, 702, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2631, 3, 6, 2547,
                                                                       663, 2553, 711, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2649, 3, 6, 2553,
                                                                       666, 2559, 720, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2667, 3, 6, 2565,
                                                                       678, 2571, 729, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2685, 3, 6, 2571,
                                                                       681, 2577, 738, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2703, 3, 6, 2577,
                                                                       684, 2583, 747, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2721, 3, 6, 2583,
                                                                       687, 2589, 756, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2739, 0, 6, 2535,
                                                                       657, 2541, 783, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2757, 0, 6, 2541,
                                                                       660, 2547, 792, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2775, 0, 6, 2547,
                                                                       663, 2553, 801, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2793, 0, 6, 2553,
                                                                       666, 2559, 810, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2811, 0, 6, 2565,
                                                                       678, 2571, 837, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2829, 0, 6, 2571,
                                                                       681, 2577, 846, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2847, 0, 6, 2577,
                                                                       684, 2583, 855, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2865, 0, 6, 2583,
                                                                       687, 2589, 864, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2883, 0, 3, 6,
                                                                       2595, 693, 2613, 2739,
                                                                       783, 2757, 97, 106, 927,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2937, 0, 3, 6,
                                                                       2613, 702, 2631, 2757,
                                                                       792, 2775, 106, 115, 954,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2991, 0, 3, 6,
                                                                       2631, 711, 2649, 2775,
                                                                       801, 2793, 115, 124, 981,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3045, 0, 3, 6,
                                                                       2667, 729, 2685, 2811,
                                                                       837, 2829, 142, 151, 1062,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3099, 0, 3, 6,
                                                                       2685, 738, 2703, 2829,
                                                                       846, 2847, 151, 160, 1089,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3153, 0, 3, 6,
                                                                       2703, 747, 2721, 2847,
                                                                       855, 2865, 160, 169, 1116,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3207, 0, 6, 2739,
                                                                       783, 2757, 187, 193, 1179,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3243, 0, 6, 2757,
                                                                       792, 2775, 193, 199, 1197,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3279, 0, 6, 2775,
                                                                       801, 2793, 199, 205, 1215,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3315, 0, 6, 2811,
                                                                       837, 2829, 217, 223, 1269,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3351, 0, 6, 2829,
                                                                       846, 2847, 223, 229, 1287,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 3387, 0, 6, 2847,
                                                                       855, 2865, 229, 235, 1305,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3423, 0, 3, 6,
                                                                       2883, 927, 2937, 3207,
                                                                       1179, 3243, 247, 265,
                                                                       1431, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3531, 0, 3, 6,
                                                                       2937, 954, 2991, 3243,
                                                                       1197, 3279, 265, 283,
                                                                       1485, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3639, 0, 3, 6,
                                                                       3045, 1062, 3099, 3315,
                                                                       1269, 3351, 319, 337,
                                                                       1647, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 3747, 0, 3, 6,
                                                                       3099, 1089, 3153, 3351,
                                                                       1287, 3387, 337, 355,
                                                                       1701, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3855, 0, 6, 3207,
                                                                       1179, 3243, 391, 401,
                                                                       1815, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3915, 0, 6, 3243,
                                                                       1197, 3279, 401, 411,
                                                                       1845, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 3975, 0, 6, 3315,
                                                                       1269, 3351, 431, 441,
                                                                       1935, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 4035, 0, 6, 3351,
                                                                       1287, 3387, 441, 451,
                                                                       1965, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 4095, 0, 3, 6,
                                                                       2883, 2937, 3423, 1431,
                                                                       3531, 3855, 1815, 3915,
                                                                       471, 501, 2175, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 4275, 0, 3, 6,
                                                                       3045, 3099, 3639, 1647,
                                                                       3747, 3975, 1935, 4035,
                                                                       561, 591, 2445, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4455, 6, 651, 654,
                                                                       2535, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4465, 6, 654, 657,
                                                                       2541, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4475, 6, 657, 660,
                                                                       2547, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4485, 6, 660, 663,
                                                                       2553, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4495, 6, 663, 666,
                                                                       2559, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4505, 6, 672, 675,
                                                                       2565, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4515, 6, 675, 678,
                                                                       2571, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4525, 6, 678, 681,
                                                                       2577, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4535, 6, 681, 684,
                                                                       2583, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4545, 6, 684, 687,
                                                                       2589, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4555, 3, 6, 4455,
                                                                       2535, 4465, 2595, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4585, 3, 6, 4465,
                                                                       2541, 4475, 2613, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4615, 3, 6, 4475,
                                                                       2547, 4485, 2631, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4645, 3, 6, 4485,
                                                                       2553, 4495, 2649, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4675, 3, 6, 4505,
                                                                       2565, 4515, 2667, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4705, 3, 6, 4515,
                                                                       2571, 4525, 2685, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4735, 3, 6, 4525,
                                                                       2577, 4535, 2703, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4765, 3, 6, 4535,
                                                                       2583, 4545, 2721, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4795, 0, 6, 4455,
                                                                       2535, 4465, 765, 774,
                                                                       2739, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4825, 0, 6, 4465,
                                                                       2541, 4475, 774, 783,
                                                                       2757, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4855, 0, 6, 4475,
                                                                       2547, 4485, 783, 792,
                                                                       2775, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4885, 0, 6, 4485,
                                                                       2553, 4495, 792, 801,
                                                                       2793, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4915, 0, 6, 4505,
                                                                       2565, 4515, 819, 828,
                                                                       2811, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4945, 0, 6, 4515,
                                                                       2571, 4525, 828, 837,
                                                                       2829, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4975, 0, 6, 4525,
                                                                       2577, 4535, 837, 846,
                                                                       2847, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5005, 0, 6, 4535,
                                                                       2583, 4545, 846, 855,
                                                                       2865, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5035, 0, 3, 6,
                                                                       4555, 2595, 4585, 4795,
                                                                       2739, 4825, 873, 900,
                                                                       2883, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5125, 0, 3, 6,
                                                                       4585, 2613, 4615, 4825,
                                                                       2757, 4855, 900, 927,
                                                                       2937, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5215, 0, 3, 6,
                                                                       4615, 2631, 4645, 4855,
                                                                       2775, 4885, 927, 954,
                                                                       2991, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5305, 0, 3, 6,
                                                                       4675, 2667, 4705, 4915,
                                                                       2811, 4945, 1008, 1035,
                                                                       3045, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5395, 0, 3, 6,
                                                                       4705, 2685, 4735, 4945,
                                                                       2829, 4975, 1035, 1062,
                                                                       3099, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5485, 0, 3, 6,
                                                                       4735, 2703, 4765, 4975,
                                                                       2847, 5005, 1062, 1089,
                                                                       3153, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5575, 0, 6, 4795,
                                                                       2739, 4825, 1143, 1161,
                                                                       3207, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5635, 0, 6, 4825,
                                                                       2757, 4855, 1161, 1179,
                                                                       3243, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5695, 0, 6, 4855,
                                                                       2775, 4885, 1179, 1197,
                                                                       3279, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5755, 0, 6, 4915,
                                                                       2811, 4945, 1233, 1251,
                                                                       3315, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5815, 0, 6, 4945,
                                                                       2829, 4975, 1251, 1269,
                                                                       3351, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 5875, 0, 6, 4975,
                                                                       2847, 5005, 1269, 1287,
                                                                       3387, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 5935, 0, 3, 6,
                                                                       5035, 2883, 5125, 5575,
                                                                       3207, 5635, 1323, 1377,
                                                                       3423, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 6115, 0, 3, 6,
                                                                       5125, 2937, 5215, 5635,
                                                                       3243, 5695, 1377, 1431,
                                                                       3531, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 6295, 0, 3, 6,
                                                                       5305, 3045, 5395, 5755,
                                                                       3315, 5815, 1539, 1593,
                                                                       3639, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 6475, 0, 3, 6,
                                                                       5395, 3099, 5485, 5815,
                                                                       3351, 5875, 1593, 1647,
                                                                       3747, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6655, 0, 6, 5575,
                                                                       3207, 5635, 1755, 1785,
                                                                       3855, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6755, 0, 6, 5635,
                                                                       3243, 5695, 1785, 1815,
                                                                       3915, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6855, 0, 6, 5755,
                                                                       3315, 5815, 1875, 1905,
                                                                       3975, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 6955, 0, 6, 5815,
                                                                       3351, 5875, 1905, 1935,
                                                                       4035, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 7055, 0, 3, 6,
                                                                       5035, 5125, 5935, 3423,
                                                                       6115, 6655, 3855, 6755,
                                                                       1995, 2085, 4095, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 7355, 0, 3, 6,
                                                                       5305, 5395, 6295, 3639,
                                                                       6475, 6855, 3975, 6955,
                                                                       2265, 2355, 4275, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 7655, 7355, 10, 10, ncols, beta);

                    simdgeo::geom_s_y(buffer, 7755, 7355, 10, 10, ncols, beta);

                    simdgeo::geom_s_z(buffer, 7855, 7355, 10, 10, ncols, beta);

                    simdgeo::geom_s_x(buffer, 7955, 7055, 10, 10, ncols, beta);

                    simdgeo::geom_s_y(buffer, 8055, 7055, 10, 10, ncols, beta);

                    simdgeo::geom_s_z(buffer, 8155, 7055, 10, 10, ncols, beta);

                    simdfunc::contract_primitives(buffer, 8255, 7655, 600, ncols);
                }
            }
        }

        simdtrf::transform_f_inner(buffer, 8855, 8255, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 8855, 7, nmax);

        simdtrf::transform_f_inner(buffer, 8855, 8355, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 49 * nvalues + n * npairs, nvalues, buffer, 8855, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 8855, 8455, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 98 * nvalues + n * npairs, nvalues, buffer, 8855, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 8855, 8555, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 147 * nvalues + n * npairs, nvalues, buffer, 8855, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 8855, 8655, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 196 * nvalues + n * npairs, nvalues, buffer, 8855, 7,
                                   nmax);

        simdtrf::transform_f_inner(buffer, 8855, 8755, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 245 * nvalues + n * npairs, nvalues, buffer, 8855, 7,
                                   nmax);
    }

    for (size_t m = 0; m < 294; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
