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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecDSD.hpp"

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

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_dsd_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_dsd_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 2033, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 150 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 2033, 1787, 216, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 9, 6, 5,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 16, 6, 5,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 23, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 26, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 29, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 32, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 3, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 3, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 3, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 53, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 56, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 59, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 62, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 65, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 68, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 71, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 74, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 77, 0, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 80, 0, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 83, 0, 6, 10, 11,
                                                                       23, 26, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 92, 0, 6, 11, 12,
                                                                       26, 29, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 101, 0, 6, 12, 13,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 110, 0, 6, 13, 14,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 119, 0, 6, 17, 18,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 128, 0, 6, 18, 19,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 137, 0, 6, 19, 20,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 146, 0, 6, 20, 21,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 155, 0, 6, 10, 11,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 161, 0, 6, 11, 12,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 167, 0, 6, 12, 13,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 173, 0, 6, 13, 14,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 179, 0, 6, 17, 18,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 185, 0, 6, 18, 19,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 191, 0, 6, 19, 20,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 197, 0, 6, 20, 21,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 203, 0, 3, 6, 53,
                                                                       56, 83, 92, 155, 161,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 221, 0, 3, 6, 56,
                                                                       59, 92, 101, 161, 167,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 239, 0, 3, 6, 59,
                                                                       62, 101, 110, 167, 173,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 257, 0, 3, 6, 68,
                                                                       71, 119, 128, 179, 185,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 275, 0, 3, 6, 71,
                                                                       74, 128, 137, 185, 191,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 293, 0, 3, 6, 74,
                                                                       77, 137, 146, 191, 197,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 311, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 314, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 317, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 320, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 323, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 326, 6, 20, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 329, 6, 21, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 332, 6, 22, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 335, 6, 12, 29,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 344, 6, 13, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 353, 6, 14, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 362, 6, 19, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 371, 6, 20, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 380, 6, 21, 50,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 389, 6, 12, 59,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 398, 6, 13, 62,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 407, 6, 14, 65,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 416, 6, 19, 74,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 425, 6, 20, 77,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 434, 6, 21, 80,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 443, 6, 29, 59,
                                                                       101, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 470, 6, 32, 62,
                                                                       110, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 497, 6, 44, 74,
                                                                       137, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 524, 6, 47, 77,
                                                                       146, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 551, 6, 59, 167,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 569, 6, 62, 173,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 587, 6, 74, 191,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 605, 6, 77, 197,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 623, 0, 6, 443,
                                                                       101, 470, 167, 239, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 677, 0, 6, 497,
                                                                       137, 524, 191, 293, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 731, 6, 10, 11,
                                                                       311, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 737, 6, 11, 12,
                                                                       314, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 743, 6, 12, 13,
                                                                       317, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 749, 6, 13, 14,
                                                                       320, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 755, 6, 17, 18,
                                                                       323, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 761, 6, 18, 19,
                                                                       326, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 767, 6, 19, 20,
                                                                       329, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 773, 6, 20, 21,
                                                                       332, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 779, 3, 6, 731,
                                                                       311, 737, 335, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 797, 3, 6, 737,
                                                                       314, 743, 344, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 815, 3, 6, 743,
                                                                       317, 749, 353, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 833, 3, 6, 755,
                                                                       323, 761, 362, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 851, 3, 6, 761,
                                                                       326, 767, 371, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 869, 3, 6, 767,
                                                                       329, 773, 380, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 887, 0, 6, 731,
                                                                       311, 737, 389, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 905, 0, 6, 737,
                                                                       314, 743, 398, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 923, 0, 6, 743,
                                                                       317, 749, 407, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 941, 0, 6, 755,
                                                                       323, 761, 416, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 959, 0, 6, 761,
                                                                       326, 767, 425, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 977, 0, 6, 767,
                                                                       329, 773, 434, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 995, 0, 3, 6, 779,
                                                                       335, 797, 887, 389, 905,
                                                                       83, 92, 443, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1049, 0, 3, 6,
                                                                       797, 344, 815, 905, 398,
                                                                       923, 92, 101, 470, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1103, 0, 3, 6,
                                                                       833, 362, 851, 941, 416,
                                                                       959, 119, 128, 497, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 1157, 0, 3, 6,
                                                                       851, 371, 869, 959, 425,
                                                                       977, 128, 137, 524, ncols,
                                                                       gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1211, 0, 6, 887,
                                                                       389, 905, 155, 161, 551,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1247, 0, 6, 905,
                                                                       398, 923, 161, 167, 569,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1283, 0, 6, 941,
                                                                       416, 959, 179, 185, 587,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 1319, 0, 6, 959,
                                                                       425, 977, 185, 191, 605,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 1355, 0, 3, 6,
                                                                       995, 443, 1049, 1211, 551,
                                                                       1247, 203, 221, 623,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 1463, 0, 3, 6,
                                                                       1103, 497, 1157, 1283,
                                                                       587, 1319, 257, 275, 677,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_s_x(buffer, 1571, 1463, 6, 6, ncols, beta);

                    simdgeo::geom_s_y(buffer, 1607, 1463, 6, 6, ncols, beta);

                    simdgeo::geom_s_z(buffer, 1643, 1463, 6, 6, ncols, beta);

                    simdgeo::geom_s_x(buffer, 1679, 1355, 6, 6, ncols, beta);

                    simdgeo::geom_s_y(buffer, 1715, 1355, 6, 6, ncols, beta);

                    simdgeo::geom_s_z(buffer, 1751, 1355, 6, 6, ncols, beta);

                    simdfunc::contract_primitives(buffer, 1787, 1571, 216, ncols);
                }
            }
        }

        simdtrf::transform_d_inner(buffer, 2003, 1787, 6, 1, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 2003, 5, nmax);

        simdtrf::transform_d_inner(buffer, 2003, 1823, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 25 * nvalues + n * npairs, nvalues, buffer, 2003, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 2003, 1859, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 50 * nvalues + n * npairs, nvalues, buffer, 2003, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 2003, 1895, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 75 * nvalues + n * npairs, nvalues, buffer, 2003, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 2003, 1931, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 100 * nvalues + n * npairs, nvalues, buffer, 2003, 5,
                                   nmax);

        simdtrf::transform_d_inner(buffer, 2003, 1967, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 125 * nvalues + n * npairs, nvalues, buffer, 2003, 5,
                                   nmax);
    }

    for (size_t m = 0; m < 150; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
