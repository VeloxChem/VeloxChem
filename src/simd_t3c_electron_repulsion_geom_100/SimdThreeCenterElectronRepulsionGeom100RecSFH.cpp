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


#include "SimdThreeCenterElectronRepulsionGeom100RecSFH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFS.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
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
compute_geom_100_sfh_three_center_electron_repulsion(double               *values,
                                                     const size_t          npairs,
                                                     const size_t          natoms,
                                                     const CBasisFunction &a_function,
                                                     const CBasisFunction &b_function,
                                                     const CBasisFunction &c_function,
                                                     const CSimdMatrix    &coordinates,
                                                     const CSimdMatrix    &c_coordinates,
                                                     CSimdMatrix          &buffer,
                                                     const double          threshold) -> void
{
    if (npairs > coordinates.number_of_columns())
    {
        errors::assertMsgCritical(
            false, std::string("compute_geom_100_sfh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 14162, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 231 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 14162, 13422, 630, dimensions);

        for (size_t i = 0; i < nprim_a; i++)
        {
            for (size_t j = 0; j < nprim_b; j++)
            {
                const auto p = a_exps[i] + b_exps[j];

                const auto alpha = a_exps[i];

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

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 9, 6, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9}, ncols, fj,
                                                        i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 19, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 43, 3, 6, 10, 11,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 49, 3, 6, 11, 12,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 55, 3, 6, 12, 13,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 61, 3, 6, 13, 14,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 67, 3, 6, 14, 15,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 73, 3, 6, 15, 16,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 79, 3, 6, 16, 17,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 85, 3, 6, 19, 22,
                                                                       43, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 95, 3, 6, 22, 25,
                                                                       49, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 105, 3, 6, 25, 28,
                                                                       55, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 115, 3, 6, 28, 31,
                                                                       61, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 125, 3, 6, 31, 34,
                                                                       67, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 135, 3, 6, 34, 37,
                                                                       73, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 145, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 148, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 151, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 154, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 157, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 160, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 163, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 166, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 169, 0, 6, 10, 11,
                                                                       19, 22, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 178, 0, 6, 11, 12,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 187, 0, 6, 12, 13,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 196, 0, 6, 13, 14,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 205, 0, 6, 14, 15,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 214, 0, 6, 15, 16,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 223, 0, 6, 16, 17,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 232, 0, 3, 6, 19,
                                                                       22, 43, 49, 169, 178,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 250, 0, 3, 6, 22,
                                                                       25, 49, 55, 178, 187,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 268, 0, 3, 6, 25,
                                                                       28, 55, 61, 187, 196,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 286, 0, 3, 6, 28,
                                                                       31, 61, 67, 196, 205,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 304, 0, 3, 6, 31,
                                                                       34, 67, 73, 205, 214,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 322, 0, 3, 6, 34,
                                                                       37, 73, 79, 214, 223,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 340, 0, 3, 6, 43,
                                                                       49, 85, 95, 232, 250,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 370, 0, 3, 6, 49,
                                                                       55, 95, 105, 250, 268,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 400, 0, 3, 6, 55,
                                                                       61, 105, 115, 268, 286,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 430, 0, 3, 6, 61,
                                                                       67, 115, 125, 286, 304,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 460, 0, 3, 6, 67,
                                                                       73, 125, 135, 304, 322,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 490, 6, 10, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 493, 6, 11, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 496, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 499, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 502, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 505, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 508, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 511, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 514, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 517, 6, 12, 25,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 526, 6, 13, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 535, 6, 14, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 544, 6, 15, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 553, 6, 16, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 562, 6, 17, 40,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 571, 6, 19, 43,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 589, 6, 22, 49,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 607, 6, 25, 55,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 625, 6, 28, 61,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 643, 6, 31, 67,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 661, 6, 34, 73,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 679, 6, 37, 79,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 697, 6, 43, 85,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 727, 6, 49, 95,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 757, 6, 55, 105,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 787, 6, 61, 115,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 817, 6, 67, 125,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 847, 6, 73, 135,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 877, 6, 10, 145,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 886, 6, 11, 148,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 895, 6, 12, 151,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 904, 6, 13, 154,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 913, 6, 14, 157,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 922, 6, 15, 160,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 931, 6, 16, 163,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 940, 6, 17, 166,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 949, 6, 19, 145,
                                                                       169, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 976, 6, 22, 148,
                                                                       178, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1003, 6, 25, 151,
                                                                       187, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1030, 6, 28, 154,
                                                                       196, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1057, 6, 31, 157,
                                                                       205, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1084, 6, 34, 160,
                                                                       214, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1111, 6, 37, 163,
                                                                       223, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1138, 3, 6, 43,
                                                                       949, 169, 976, 232, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1192, 3, 6, 49,
                                                                       976, 178, 1003, 250,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1246, 3, 6, 55,
                                                                       1003, 187, 1030, 268,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1300, 3, 6, 61,
                                                                       1030, 196, 1057, 286,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1354, 3, 6, 67,
                                                                       1057, 205, 1084, 304,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1408, 3, 6, 73,
                                                                       1084, 214, 1111, 322,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1462, 3, 6, 85,
                                                                       1138, 232, 1192, 340,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1552, 3, 6, 95,
                                                                       1192, 250, 1246, 370,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1642, 3, 6, 105,
                                                                       1246, 268, 1300, 400,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1732, 3, 6, 115,
                                                                       1300, 286, 1354, 430,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 1822, 3, 6, 125,
                                                                       1354, 304, 1408, 460,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1912, 6, 10, 11,
                                                                       496, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1918, 6, 11, 12,
                                                                       499, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1924, 6, 12, 13,
                                                                       502, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1930, 6, 13, 14,
                                                                       505, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1936, 6, 14, 15,
                                                                       508, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1942, 6, 15, 16,
                                                                       511, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1948, 6, 16, 17,
                                                                       514, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1954, 3, 6, 1912,
                                                                       496, 1918, 517, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1972, 3, 6, 1918,
                                                                       499, 1924, 526, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 1990, 3, 6, 1924,
                                                                       502, 1930, 535, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2008, 3, 6, 1930,
                                                                       505, 1936, 544, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2026, 3, 6, 1936,
                                                                       508, 1942, 553, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2044, 3, 6, 1942,
                                                                       511, 1948, 562, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2062, 3, 6, 1954,
                                                                       517, 1972, 43, 49, 607,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2098, 3, 6, 1972,
                                                                       526, 1990, 49, 55, 625,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2134, 3, 6, 1990,
                                                                       535, 2008, 55, 61, 643,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2170, 3, 6, 2008,
                                                                       544, 2026, 61, 67, 661,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2206, 3, 6, 2026,
                                                                       553, 2044, 67, 73, 679,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2242, 3, 6, 2062,
                                                                       607, 2098, 85, 95, 757,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2302, 3, 6, 2098,
                                                                       625, 2134, 95, 105, 787,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2362, 3, 6, 2134,
                                                                       643, 2170, 105, 115, 817,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 2422, 3, 6, 2170,
                                                                       661, 2206, 115, 125, 847,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2482, 0, 6, 1912,
                                                                       496, 1918, 895, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2500, 0, 6, 1918,
                                                                       499, 1924, 904, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2518, 0, 6, 1924,
                                                                       502, 1930, 913, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2536, 0, 6, 1930,
                                                                       505, 1936, 922, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2554, 0, 6, 1936,
                                                                       508, 1942, 931, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2572, 0, 6, 1942,
                                                                       511, 1948, 940, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2590, 0, 3, 6,
                                                                       1954, 517, 1972, 2482,
                                                                       895, 2500, 169, 178, 1003,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2644, 0, 3, 6,
                                                                       1972, 526, 1990, 2500,
                                                                       904, 2518, 178, 187, 1030,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2698, 0, 3, 6,
                                                                       1990, 535, 2008, 2518,
                                                                       913, 2536, 187, 196, 1057,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2752, 0, 3, 6,
                                                                       2008, 544, 2026, 2536,
                                                                       922, 2554, 196, 205, 1084,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 2806, 0, 3, 6,
                                                                       2026, 553, 2044, 2554,
                                                                       931, 2572, 205, 214, 1111,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 2860, 0, 3, 6,
                                                                       2062, 607, 2098, 2590,
                                                                       1003, 2644, 232, 250,
                                                                       1246, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 2968, 0, 3, 6,
                                                                       2098, 625, 2134, 2644,
                                                                       1030, 2698, 250, 268,
                                                                       1300, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3076, 0, 3, 6,
                                                                       2134, 643, 2170, 2698,
                                                                       1057, 2752, 268, 286,
                                                                       1354, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3184, 0, 3, 6,
                                                                       2170, 661, 2206, 2752,
                                                                       1084, 2806, 286, 304,
                                                                       1408, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 3292, 0, 3, 6,
                                                                       2242, 757, 2302, 2590,
                                                                       2644, 2860, 1246, 2968,
                                                                       340, 370, 1642, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 3472, 0, 3, 6,
                                                                       2302, 787, 2362, 2644,
                                                                       2698, 2968, 1300, 3076,
                                                                       370, 400, 1732, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 3652, 0, 3, 6,
                                                                       2362, 817, 2422, 2698,
                                                                       2752, 3076, 1354, 3184,
                                                                       400, 430, 1822, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3832, 6, 490, 493,
                                                                       1912, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3842, 6, 493, 496,
                                                                       1918, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3852, 6, 496, 499,
                                                                       1924, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3862, 6, 499, 502,
                                                                       1930, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3872, 6, 502, 505,
                                                                       1936, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3882, 6, 505, 508,
                                                                       1942, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 3892, 6, 508, 511,
                                                                       1948, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3902, 3, 6, 3832,
                                                                       1912, 3842, 1954, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3932, 3, 6, 3842,
                                                                       1918, 3852, 1972, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3962, 3, 6, 3852,
                                                                       1924, 3862, 1990, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 3992, 3, 6, 3862,
                                                                       1930, 3872, 2008, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4022, 3, 6, 3872,
                                                                       1936, 3882, 2026, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4052, 3, 6, 3882,
                                                                       1942, 3892, 2044, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4082, 3, 6, 3902,
                                                                       1954, 3932, 571, 589,
                                                                       2062, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4142, 3, 6, 3932,
                                                                       1972, 3962, 589, 607,
                                                                       2098, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4202, 3, 6, 3962,
                                                                       1990, 3992, 607, 625,
                                                                       2134, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4262, 3, 6, 3992,
                                                                       2008, 4022, 625, 643,
                                                                       2170, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 4322, 3, 6, 4022,
                                                                       2026, 4052, 643, 661,
                                                                       2206, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4382, 3, 6, 4082,
                                                                       2062, 4142, 697, 727,
                                                                       2242, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4482, 3, 6, 4142,
                                                                       2098, 4202, 727, 757,
                                                                       2302, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4582, 3, 6, 4202,
                                                                       2134, 4262, 757, 787,
                                                                       2362, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 4682, 3, 6, 4262,
                                                                       2170, 4322, 787, 817,
                                                                       2422, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4782, 0, 6, 3832,
                                                                       1912, 3842, 877, 886,
                                                                       2482, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4812, 0, 6, 3842,
                                                                       1918, 3852, 886, 895,
                                                                       2500, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4842, 0, 6, 3852,
                                                                       1924, 3862, 895, 904,
                                                                       2518, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4872, 0, 6, 3862,
                                                                       1930, 3872, 904, 913,
                                                                       2536, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4902, 0, 6, 3872,
                                                                       1936, 3882, 913, 922,
                                                                       2554, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 4932, 0, 6, 3882,
                                                                       1942, 3892, 922, 931,
                                                                       2572, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 4962, 0, 3, 6,
                                                                       3902, 1954, 3932, 4782,
                                                                       2482, 4812, 949, 976,
                                                                       2590, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5052, 0, 3, 6,
                                                                       3932, 1972, 3962, 4812,
                                                                       2500, 4842, 976, 1003,
                                                                       2644, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5142, 0, 3, 6,
                                                                       3962, 1990, 3992, 4842,
                                                                       2518, 4872, 1003, 1030,
                                                                       2698, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5232, 0, 3, 6,
                                                                       3992, 2008, 4022, 4872,
                                                                       2536, 4902, 1030, 1057,
                                                                       2752, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5322, 0, 3, 6,
                                                                       4022, 2026, 4052, 4902,
                                                                       2554, 4932, 1057, 1084,
                                                                       2806, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 5412, 0, 3, 6,
                                                                       4082, 2062, 4142, 4962,
                                                                       2590, 5052, 1138, 1192,
                                                                       2860, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 5592, 0, 3, 6,
                                                                       4142, 2098, 4202, 5052,
                                                                       2644, 5142, 1192, 1246,
                                                                       2968, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 5772, 0, 3, 6,
                                                                       4202, 2134, 4262, 5142,
                                                                       2698, 5232, 1246, 1300,
                                                                       3076, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 5952, 0, 3, 6,
                                                                       4262, 2170, 4322, 5232,
                                                                       2752, 5322, 1300, 1354,
                                                                       3184, ncols, gamma, p,
                                                                       q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 6132, 0, 3, 6,
                                                                       4382, 2242, 4482, 4962,
                                                                       5052, 5412, 2860, 5592,
                                                                       1462, 1552, 3292, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 6432, 0, 3, 6,
                                                                       4482, 2302, 4582, 5052,
                                                                       5142, 5592, 2968, 5772,
                                                                       1552, 1642, 3472, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 6732, 0, 3, 6,
                                                                       4582, 2362, 4682, 5142,
                                                                       5232, 5772, 3076, 5952,
                                                                       1642, 1732, 3652, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7032, 6, 1912,
                                                                       1918, 3852, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7047, 6, 1918,
                                                                       1924, 3862, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7062, 6, 1924,
                                                                       1930, 3872, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7077, 6, 1930,
                                                                       1936, 3882, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7092, 6, 1936,
                                                                       1942, 3892, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7107, 3, 6, 7032,
                                                                       3852, 7047, 1954, 1972,
                                                                       3962, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7152, 3, 6, 7047,
                                                                       3862, 7062, 1972, 1990,
                                                                       3992, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7197, 3, 6, 7062,
                                                                       3872, 7077, 1990, 2008,
                                                                       4022, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7242, 3, 6, 7077,
                                                                       3882, 7092, 2008, 2026,
                                                                       4052, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7287, 3, 6, 7107,
                                                                       3962, 7152, 2062, 2098,
                                                                       4202, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7377, 3, 6, 7152,
                                                                       3992, 7197, 2098, 2134,
                                                                       4262, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 7467, 3, 6, 7197,
                                                                       4022, 7242, 2134, 2170,
                                                                       4322, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 7557, 3, 6, 7287,
                                                                       4202, 7377, 2242, 2302,
                                                                       4582, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 7707, 3, 6, 7377,
                                                                       4262, 7467, 2302, 2362,
                                                                       4682, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7857, 0, 6, 7032,
                                                                       3852, 7047, 2482, 2500,
                                                                       4842, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7902, 0, 6, 7047,
                                                                       3862, 7062, 2500, 2518,
                                                                       4872, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7947, 0, 6, 7062,
                                                                       3872, 7077, 2518, 2536,
                                                                       4902, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 7992, 0, 6, 7077,
                                                                       3882, 7092, 2536, 2554,
                                                                       4932, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 8037, 0, 3, 6,
                                                                       7107, 3962, 7152, 7857,
                                                                       4842, 7902, 2590, 2644,
                                                                       5142, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 8172, 0, 3, 6,
                                                                       7152, 3992, 7197, 7902,
                                                                       4872, 7947, 2644, 2698,
                                                                       5232, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 8307, 0, 3, 6,
                                                                       7197, 4022, 7242, 7947,
                                                                       4902, 7992, 2698, 2752,
                                                                       5322, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 8442, 0, 3, 6,
                                                                       7287, 4202, 7377, 8037,
                                                                       5142, 8172, 2860, 2968,
                                                                       5772, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 8712, 0, 3, 6,
                                                                       7377, 4262, 7467, 8172,
                                                                       5232, 8307, 2968, 3076,
                                                                       5952, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 8982, 0, 3, 6,
                                                                       7557, 4582, 7707, 8037,
                                                                       8172, 8442, 5772, 8712,
                                                                       3292, 3472, 6732, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9432, 6, 3832,
                                                                       3842, 7032, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9453, 6, 3842,
                                                                       3852, 7047, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9474, 6, 3852,
                                                                       3862, 7062, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9495, 6, 3862,
                                                                       3872, 7077, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 9516, 6, 3872,
                                                                       3882, 7092, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9537, 3, 6, 9432,
                                                                       7032, 9453, 3902, 3932,
                                                                       7107, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9600, 3, 6, 9453,
                                                                       7047, 9474, 3932, 3962,
                                                                       7152, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9663, 3, 6, 9474,
                                                                       7062, 9495, 3962, 3992,
                                                                       7197, ncols, gamma, p,
                                                                       q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 9726, 3, 6, 9495,
                                                                       7077, 9516, 3992, 4022,
                                                                       7242, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 9789, 3, 6, 9537,
                                                                       7107, 9600, 4082, 4142,
                                                                       7287, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 9915, 3, 6, 9600,
                                                                       7152, 9663, 4142, 4202,
                                                                       7377, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 10041, 3, 6, 9663,
                                                                       7197, 9726, 4202, 4262,
                                                                       7467, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 10167, 3, 6, 9789,
                                                                       7287, 9915, 4382, 4482,
                                                                       7557, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 10377, 3, 6, 9915,
                                                                       7377, 10041, 4482, 4582,
                                                                       7707, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 10587, 0, 6, 9432,
                                                                       7032, 9453, 4782, 4812,
                                                                       7857, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 10650, 0, 6, 9453,
                                                                       7047, 9474, 4812, 4842,
                                                                       7902, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 10713, 0, 6, 9474,
                                                                       7062, 9495, 4842, 4872,
                                                                       7947, ncols, gamma, p,
                                                                       q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 10776, 0, 6, 9495,
                                                                       7077, 9516, 4872, 4902,
                                                                       7992, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 10839, 0, 3, 6,
                                                                       9537, 7107, 9600, 10587,
                                                                       7857, 10650, 4962, 5052,
                                                                       8037, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 11028, 0, 3, 6,
                                                                       9600, 7152, 9663, 10650,
                                                                       7902, 10713, 5052, 5142,
                                                                       8172, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 11217, 0, 3, 6,
                                                                       9663, 7197, 9726, 10713,
                                                                       7947, 10776, 5142, 5232,
                                                                       8307, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 11406, 0, 3, 6,
                                                                       9789, 7287, 9915, 10839,
                                                                       8037, 11028, 5412, 5592,
                                                                       8442, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 11784, 0, 3, 6,
                                                                       9915, 7377, 10041, 11028,
                                                                       8172, 11217, 5592, 5772,
                                                                       8712, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 12162, 0, 3, 6,
                                                                       10167, 7557, 10377, 10839,
                                                                       11028, 11406, 8442, 11784,
                                                                       6132, 6432, 8982, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 12792, 12162, 1, 210, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 13002, 12162, 1, 210, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 13212, 12162, 1, 210, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 13422, 12792, 630, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 14052, 13422, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 14052, 11, nmax);

        simdtrf::transform_h_inner(buffer, 14052, 13632, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 77 * nvalues + n * npairs, nvalues, buffer, 14052,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 14052, 13842, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 154 * nvalues + n * npairs, nvalues, buffer, 14052,
                                   11, nmax);
    }

    for (size_t m = 0; m < 231; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
