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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecSFG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFS.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
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
compute_rs_geom_100_sfg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_sfg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 16609, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 16609, 15619, 900, dimensions);

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

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 77, 3, 6, 10, 11,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 83, 3, 6, 11, 12,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 89, 3, 6, 12, 13,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 95, 3, 6, 13, 14,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 101, 3, 6, 14, 15,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 107, 3, 6, 15, 16,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 113, 3, 6, 16, 17,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 119, 3, 6, 20, 21,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 125, 3, 6, 21, 22,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 131, 3, 6, 22, 23,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 137, 3, 6, 23, 24,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 143, 3, 6, 24, 25,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 149, 3, 6, 25, 26,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 155, 3, 6, 26, 27,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 161, 3, 6, 29, 32,
                                                                       77, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 171, 3, 6, 32, 35,
                                                                       83, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 181, 3, 6, 35, 38,
                                                                       89, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 191, 3, 6, 38, 41,
                                                                       95, 101, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 201, 3, 6, 41, 44,
                                                                       101, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 211, 3, 6, 44, 47,
                                                                       107, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 221, 3, 6, 53, 56,
                                                                       119, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 231, 3, 6, 56, 59,
                                                                       125, 131, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 241, 3, 6, 59, 62,
                                                                       131, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 251, 3, 6, 62, 65,
                                                                       137, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 261, 3, 6, 65, 68,
                                                                       143, 149, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 271, 3, 6, 68, 71,
                                                                       149, 155, ncols, gamma, p,
                                                                       q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 281, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 284, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 287, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 290, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 293, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 296, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 299, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 302, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 305, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 308, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 311, 0, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 314, 0, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 317, 0, 6, 10, 11,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 326, 0, 6, 11, 12,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 335, 0, 6, 12, 13,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 344, 0, 6, 13, 14,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 353, 0, 6, 14, 15,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 362, 0, 6, 15, 16,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 371, 0, 6, 16, 17,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 380, 0, 6, 20, 21,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 389, 0, 6, 21, 22,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 398, 0, 6, 22, 23,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 407, 0, 6, 23, 24,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 416, 0, 6, 24, 25,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 425, 0, 6, 25, 26,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 434, 0, 6, 26, 27,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 443, 0, 3, 6, 29,
                                                                       32, 77, 83, 317, 326,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 461, 0, 3, 6, 32,
                                                                       35, 83, 89, 326, 335,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 479, 0, 3, 6, 35,
                                                                       38, 89, 95, 335, 344,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 497, 0, 3, 6, 38,
                                                                       41, 95, 101, 344, 353,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 515, 0, 3, 6, 41,
                                                                       44, 101, 107, 353, 362,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 533, 0, 3, 6, 44,
                                                                       47, 107, 113, 362, 371,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 551, 0, 3, 6, 53,
                                                                       56, 119, 125, 380, 389,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 569, 0, 3, 6, 56,
                                                                       59, 125, 131, 389, 398,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 587, 0, 3, 6, 59,
                                                                       62, 131, 137, 398, 407,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 605, 0, 3, 6, 62,
                                                                       65, 137, 143, 407, 416,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 623, 0, 3, 6, 65,
                                                                       68, 143, 149, 416, 425,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 641, 0, 3, 6, 68,
                                                                       71, 149, 155, 425, 434,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 659, 0, 3, 6, 77,
                                                                       83, 161, 171, 443, 461,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 689, 0, 3, 6, 83,
                                                                       89, 171, 181, 461, 479,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 719, 0, 3, 6, 89,
                                                                       95, 181, 191, 479, 497,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 749, 0, 3, 6, 95,
                                                                       101, 191, 201, 497, 515,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 779, 0, 3, 6, 101,
                                                                       107, 201, 211, 515, 533,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 809, 0, 3, 6, 119,
                                                                       125, 221, 231, 551, 569,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 839, 0, 3, 6, 125,
                                                                       131, 231, 241, 569, 587,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 869, 0, 3, 6, 131,
                                                                       137, 241, 251, 587, 605,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 899, 0, 3, 6, 137,
                                                                       143, 251, 261, 605, 623,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 929, 0, 3, 6, 143,
                                                                       149, 261, 271, 623, 641,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 959, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 962, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 965, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 968, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 971, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 974, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 977, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 980, 6, 22, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 983, 6, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 986, 6, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 989, 6, 25, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 992, 6, 26, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 995, 6, 27, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 998, 6, 28, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1001, 6, 12, 35,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1010, 6, 13, 38,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1019, 6, 14, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1028, 6, 15, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1037, 6, 16, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1046, 6, 17, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1055, 6, 22, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1064, 6, 23, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1073, 6, 24, 65,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1082, 6, 25, 68,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1091, 6, 26, 71,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1100, 6, 27, 74,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1109, 6, 35, 89,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1127, 6, 38, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1145, 6, 41, 101,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1163, 6, 44, 107,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1181, 6, 47, 113,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1199, 6, 59, 131,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1217, 6, 62, 137,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1235, 6, 65, 143,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1253, 6, 68, 149,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1271, 6, 71, 155,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1289, 6, 89, 181,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1319, 6, 95, 191,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1349, 6, 101, 201,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1379, 6, 107, 211,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1409, 6, 131, 241,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1439, 6, 137, 251,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1469, 6, 143, 261,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1499, 6, 149, 271,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1529, 6, 12, 281,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1538, 6, 13, 284,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1547, 6, 14, 287,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1556, 6, 15, 290,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1565, 6, 16, 293,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1574, 6, 17, 296,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1583, 6, 22, 299,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1592, 6, 23, 302,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1601, 6, 24, 305,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1610, 6, 25, 308,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1619, 6, 26, 311,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1628, 6, 27, 314,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1637, 6, 35, 281,
                                                                       335, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1664, 6, 38, 284,
                                                                       344, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1691, 6, 41, 287,
                                                                       353, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1718, 6, 44, 290,
                                                                       362, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1745, 6, 47, 293,
                                                                       371, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1772, 6, 59, 299,
                                                                       398, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1799, 6, 62, 302,
                                                                       407, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1826, 6, 65, 305,
                                                                       416, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1853, 6, 68, 308,
                                                                       425, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1880, 6, 71, 311,
                                                                       434, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1907, 3, 6, 89,
                                                                       1637, 335, 1664, 479,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1961, 3, 6, 95,
                                                                       1664, 344, 1691, 497,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2015, 3, 6, 101,
                                                                       1691, 353, 1718, 515,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2069, 3, 6, 107,
                                                                       1718, 362, 1745, 533,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2123, 3, 6, 131,
                                                                       1772, 398, 1799, 587,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2177, 3, 6, 137,
                                                                       1799, 407, 1826, 605,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2231, 3, 6, 143,
                                                                       1826, 416, 1853, 623,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2285, 3, 6, 149,
                                                                       1853, 425, 1880, 641,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2339, 3, 6, 181,
                                                                       1907, 479, 1961, 719,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2429, 3, 6, 191,
                                                                       1961, 497, 2015, 749,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2519, 3, 6, 201,
                                                                       2015, 515, 2069, 779,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2609, 3, 6, 241,
                                                                       2123, 587, 2177, 869,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2699, 3, 6, 251,
                                                                       2177, 605, 2231, 899,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2789, 3, 6, 261,
                                                                       2231, 623, 2285, 929,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2879, 6, 10, 11,
                                                                       959, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2885, 6, 11, 12,
                                                                       962, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2891, 6, 12, 13,
                                                                       965, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2897, 6, 13, 14,
                                                                       968, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2903, 6, 14, 15,
                                                                       971, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2909, 6, 15, 16,
                                                                       974, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2915, 6, 16, 17,
                                                                       977, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2921, 6, 20, 21,
                                                                       980, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2927, 6, 21, 22,
                                                                       983, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2933, 6, 22, 23,
                                                                       986, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2939, 6, 23, 24,
                                                                       989, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2945, 6, 24, 25,
                                                                       992, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2951, 6, 25, 26,
                                                                       995, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2957, 6, 26, 27,
                                                                       998, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2963, 3, 6, 2879,
                                                                       959, 2885, 1001, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2981, 3, 6, 2885,
                                                                       962, 2891, 1010, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2999, 3, 6, 2891,
                                                                       965, 2897, 1019, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3017, 3, 6, 2897,
                                                                       968, 2903, 1028, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3035, 3, 6, 2903,
                                                                       971, 2909, 1037, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3053, 3, 6, 2909,
                                                                       974, 2915, 1046, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3071, 3, 6, 2921,
                                                                       980, 2927, 1055, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3089, 3, 6, 2927,
                                                                       983, 2933, 1064, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3107, 3, 6, 2933,
                                                                       986, 2939, 1073, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3125, 3, 6, 2939,
                                                                       989, 2945, 1082, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3143, 3, 6, 2945,
                                                                       992, 2951, 1091, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3161, 3, 6, 2951,
                                                                       995, 2957, 1100, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3179, 3, 6, 2963,
                                                                       1001, 2981, 77, 83, 1109,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3215, 3, 6, 2981,
                                                                       1010, 2999, 83, 89, 1127,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3251, 3, 6, 2999,
                                                                       1019, 3017, 89, 95, 1145,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3287, 3, 6, 3017,
                                                                       1028, 3035, 95, 101, 1163,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3323, 3, 6, 3035,
                                                                       1037, 3053, 101, 107,
                                                                       1181, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3359, 3, 6, 3071,
                                                                       1055, 3089, 119, 125,
                                                                       1199, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3395, 3, 6, 3089,
                                                                       1064, 3107, 125, 131,
                                                                       1217, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3431, 3, 6, 3107,
                                                                       1073, 3125, 131, 137,
                                                                       1235, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3467, 3, 6, 3125,
                                                                       1082, 3143, 137, 143,
                                                                       1253, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 3503, 3, 6, 3143,
                                                                       1091, 3161, 143, 149,
                                                                       1271, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3539, 3, 6, 3179,
                                                                       1109, 3215, 161, 171,
                                                                       1289, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3599, 3, 6, 3215,
                                                                       1127, 3251, 171, 181,
                                                                       1319, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3659, 3, 6, 3251,
                                                                       1145, 3287, 181, 191,
                                                                       1349, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3719, 3, 6, 3287,
                                                                       1163, 3323, 191, 201,
                                                                       1379, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3779, 3, 6, 3359,
                                                                       1199, 3395, 221, 231,
                                                                       1409, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3839, 3, 6, 3395,
                                                                       1217, 3431, 231, 241,
                                                                       1439, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3899, 3, 6, 3431,
                                                                       1235, 3467, 241, 251,
                                                                       1469, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 3959, 3, 6, 3467,
                                                                       1253, 3503, 251, 261,
                                                                       1499, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4019, 0, 6, 2879,
                                                                       959, 2885, 1529, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4037, 0, 6, 2885,
                                                                       962, 2891, 1538, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4055, 0, 6, 2891,
                                                                       965, 2897, 1547, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4073, 0, 6, 2897,
                                                                       968, 2903, 1556, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4091, 0, 6, 2903,
                                                                       971, 2909, 1565, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4109, 0, 6, 2909,
                                                                       974, 2915, 1574, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4127, 0, 6, 2921,
                                                                       980, 2927, 1583, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4145, 0, 6, 2927,
                                                                       983, 2933, 1592, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4163, 0, 6, 2933,
                                                                       986, 2939, 1601, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4181, 0, 6, 2939,
                                                                       989, 2945, 1610, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4199, 0, 6, 2945,
                                                                       992, 2951, 1619, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4217, 0, 6, 2951,
                                                                       995, 2957, 1628, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4235, 0, 3, 6,
                                                                       2963, 1001, 2981, 4019,
                                                                       1529, 4037, 317, 326,
                                                                       1637, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4289, 0, 3, 6,
                                                                       2981, 1010, 2999, 4037,
                                                                       1538, 4055, 326, 335,
                                                                       1664, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4343, 0, 3, 6,
                                                                       2999, 1019, 3017, 4055,
                                                                       1547, 4073, 335, 344,
                                                                       1691, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4397, 0, 3, 6,
                                                                       3017, 1028, 3035, 4073,
                                                                       1556, 4091, 344, 353,
                                                                       1718, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4451, 0, 3, 6,
                                                                       3035, 1037, 3053, 4091,
                                                                       1565, 4109, 353, 362,
                                                                       1745, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4505, 0, 3, 6,
                                                                       3071, 1055, 3089, 4127,
                                                                       1583, 4145, 380, 389,
                                                                       1772, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4559, 0, 3, 6,
                                                                       3089, 1064, 3107, 4145,
                                                                       1592, 4163, 389, 398,
                                                                       1799, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4613, 0, 3, 6,
                                                                       3107, 1073, 3125, 4163,
                                                                       1601, 4181, 398, 407,
                                                                       1826, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4667, 0, 3, 6,
                                                                       3125, 1082, 3143, 4181,
                                                                       1610, 4199, 407, 416,
                                                                       1853, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4721, 0, 3, 6,
                                                                       3143, 1091, 3161, 4199,
                                                                       1619, 4217, 416, 425,
                                                                       1880, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4775, 0, 3, 6,
                                                                       3179, 1109, 3215, 4235,
                                                                       1637, 4289, 443, 461,
                                                                       1907, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4883, 0, 3, 6,
                                                                       3215, 1127, 3251, 4289,
                                                                       1664, 4343, 461, 479,
                                                                       1961, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4991, 0, 3, 6,
                                                                       3251, 1145, 3287, 4343,
                                                                       1691, 4397, 479, 497,
                                                                       2015, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5099, 0, 3, 6,
                                                                       3287, 1163, 3323, 4397,
                                                                       1718, 4451, 497, 515,
                                                                       2069, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5207, 0, 3, 6,
                                                                       3359, 1199, 3395, 4505,
                                                                       1772, 4559, 551, 569,
                                                                       2123, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5315, 0, 3, 6,
                                                                       3395, 1217, 3431, 4559,
                                                                       1799, 4613, 569, 587,
                                                                       2177, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5423, 0, 3, 6,
                                                                       3431, 1235, 3467, 4613,
                                                                       1826, 4667, 587, 605,
                                                                       2231, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5531, 0, 3, 6,
                                                                       3467, 1253, 3503, 4667,
                                                                       1853, 4721, 605, 623,
                                                                       2285, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 5639, 0, 3, 6,
                                                                       3539, 1289, 3599, 4235,
                                                                       4289, 4775, 1907, 4883,
                                                                       659, 689, 2339, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 5819, 0, 3, 6,
                                                                       3599, 1319, 3659, 4289,
                                                                       4343, 4883, 1961, 4991,
                                                                       689, 719, 2429, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 5999, 0, 3, 6,
                                                                       3659, 1349, 3719, 4343,
                                                                       4397, 4991, 2015, 5099,
                                                                       719, 749, 2519, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 6179, 0, 3, 6,
                                                                       3779, 1409, 3839, 4505,
                                                                       4559, 5207, 2123, 5315,
                                                                       809, 839, 2609, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 6359, 0, 3, 6,
                                                                       3839, 1439, 3899, 4559,
                                                                       4613, 5315, 2177, 5423,
                                                                       839, 869, 2699, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 6539, 0, 3, 6,
                                                                       3899, 1469, 3959, 4613,
                                                                       4667, 5423, 2231, 5531,
                                                                       869, 899, 2789, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6719, 6, 959, 962,
                                                                       2891, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6729, 6, 962, 965,
                                                                       2897, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6739, 6, 965, 968,
                                                                       2903, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6749, 6, 968, 971,
                                                                       2909, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6759, 6, 971, 974,
                                                                       2915, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6769, 6, 980, 983,
                                                                       2933, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6779, 6, 983, 986,
                                                                       2939, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6789, 6, 986, 989,
                                                                       2945, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6799, 6, 989, 992,
                                                                       2951, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 6809, 6, 992, 995,
                                                                       2957, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6819, 3, 6, 6719,
                                                                       2891, 6729, 2999, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6849, 3, 6, 6729,
                                                                       2897, 6739, 3017, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6879, 3, 6, 6739,
                                                                       2903, 6749, 3035, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6909, 3, 6, 6749,
                                                                       2909, 6759, 3053, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6939, 3, 6, 6769,
                                                                       2933, 6779, 3107, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6969, 3, 6, 6779,
                                                                       2939, 6789, 3125, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 6999, 3, 6, 6789,
                                                                       2945, 6799, 3143, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 7029, 3, 6, 6799,
                                                                       2951, 6809, 3161, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7059, 3, 6, 6819,
                                                                       2999, 6849, 1109, 1127,
                                                                       3251, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7119, 3, 6, 6849,
                                                                       3017, 6879, 1127, 1145,
                                                                       3287, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7179, 3, 6, 6879,
                                                                       3035, 6909, 1145, 1163,
                                                                       3323, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7239, 3, 6, 6939,
                                                                       3107, 6969, 1199, 1217,
                                                                       3431, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7299, 3, 6, 6969,
                                                                       3125, 6999, 1217, 1235,
                                                                       3467, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 7359, 3, 6, 6999,
                                                                       3143, 7029, 1235, 1253,
                                                                       3503, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7419, 3, 6, 7059,
                                                                       3251, 7119, 1289, 1319,
                                                                       3659, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7519, 3, 6, 7119,
                                                                       3287, 7179, 1319, 1349,
                                                                       3719, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7619, 3, 6, 7239,
                                                                       3431, 7299, 1409, 1439,
                                                                       3899, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 7719, 3, 6, 7299,
                                                                       3467, 7359, 1439, 1469,
                                                                       3959, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7819, 0, 6, 6719,
                                                                       2891, 6729, 1529, 1538,
                                                                       4055, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7849, 0, 6, 6729,
                                                                       2897, 6739, 1538, 1547,
                                                                       4073, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7879, 0, 6, 6739,
                                                                       2903, 6749, 1547, 1556,
                                                                       4091, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7909, 0, 6, 6749,
                                                                       2909, 6759, 1556, 1565,
                                                                       4109, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7939, 0, 6, 6769,
                                                                       2933, 6779, 1583, 1592,
                                                                       4163, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7969, 0, 6, 6779,
                                                                       2939, 6789, 1592, 1601,
                                                                       4181, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 7999, 0, 6, 6789,
                                                                       2945, 6799, 1601, 1610,
                                                                       4199, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 8029, 0, 6, 6799,
                                                                       2951, 6809, 1610, 1619,
                                                                       4217, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 8059, 0, 3, 6,
                                                                       6819, 2999, 6849, 7819,
                                                                       4055, 7849, 1637, 1664,
                                                                       4343, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 8149, 0, 3, 6,
                                                                       6849, 3017, 6879, 7849,
                                                                       4073, 7879, 1664, 1691,
                                                                       4397, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 8239, 0, 3, 6,
                                                                       6879, 3035, 6909, 7879,
                                                                       4091, 7909, 1691, 1718,
                                                                       4451, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 8329, 0, 3, 6,
                                                                       6939, 3107, 6969, 7939,
                                                                       4163, 7969, 1772, 1799,
                                                                       4613, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 8419, 0, 3, 6,
                                                                       6969, 3125, 6999, 7969,
                                                                       4181, 7999, 1799, 1826,
                                                                       4667, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 8509, 0, 3, 6,
                                                                       6999, 3143, 7029, 7999,
                                                                       4199, 8029, 1826, 1853,
                                                                       4721, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 8599, 0, 3, 6,
                                                                       7059, 3251, 7119, 8059,
                                                                       4343, 8149, 1907, 1961,
                                                                       4991, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 8779, 0, 3, 6,
                                                                       7119, 3287, 7179, 8149,
                                                                       4397, 8239, 1961, 2015,
                                                                       5099, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 8959, 0, 3, 6,
                                                                       7239, 3431, 7299, 8329,
                                                                       4613, 8419, 2123, 2177,
                                                                       5423, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 9139, 0, 3, 6,
                                                                       7299, 3467, 7359, 8419,
                                                                       4667, 8509, 2177, 2231,
                                                                       5531, ncols, gamma, p,
                                                                       q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 9319, 0, 3, 6,
                                                                       7419, 3659, 7519, 8059,
                                                                       8149, 8599, 4991, 8779,
                                                                       2339, 2429, 5999, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 9619, 0, 3, 6,
                                                                       7619, 3899, 7719, 8329,
                                                                       8419, 8959, 5423, 9139,
                                                                       2609, 2699, 6539, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9919, 6, 2879,
                                                                       2885, 6719, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9934, 6, 2885,
                                                                       2891, 6729, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9949, 6, 2891,
                                                                       2897, 6739, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9964, 6, 2897,
                                                                       2903, 6749, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9979, 6, 2903,
                                                                       2909, 6759, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 9994, 6, 2921,
                                                                       2927, 6769, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10009, 6, 2927,
                                                                       2933, 6779, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10024, 6, 2933,
                                                                       2939, 6789, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10039, 6, 2939,
                                                                       2945, 6799, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 10054, 6, 2945,
                                                                       2951, 6809, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10069, 3, 6, 9919,
                                                                       6719, 9934, 2963, 2981,
                                                                       6819, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10114, 3, 6, 9934,
                                                                       6729, 9949, 2981, 2999,
                                                                       6849, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10159, 3, 6, 9949,
                                                                       6739, 9964, 2999, 3017,
                                                                       6879, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10204, 3, 6, 9964,
                                                                       6749, 9979, 3017, 3035,
                                                                       6909, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10249, 3, 6, 9994,
                                                                       6769, 10009, 3071, 3089,
                                                                       6939, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10294, 3, 6,
                                                                       10009, 6779, 10024, 3089,
                                                                       3107, 6969, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10339, 3, 6,
                                                                       10024, 6789, 10039, 3107,
                                                                       3125, 6999, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 10384, 3, 6,
                                                                       10039, 6799, 10054, 3125,
                                                                       3143, 7029, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 10429, 3, 6,
                                                                       10069, 6819, 10114, 3179,
                                                                       3215, 7059, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 10519, 3, 6,
                                                                       10114, 6849, 10159, 3215,
                                                                       3251, 7119, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 10609, 3, 6,
                                                                       10159, 6879, 10204, 3251,
                                                                       3287, 7179, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 10699, 3, 6,
                                                                       10249, 6939, 10294, 3359,
                                                                       3395, 7239, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 10789, 3, 6,
                                                                       10294, 6969, 10339, 3395,
                                                                       3431, 7299, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 10879, 3, 6,
                                                                       10339, 6999, 10384, 3431,
                                                                       3467, 7359, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 10969, 3, 6,
                                                                       10429, 7059, 10519, 3539,
                                                                       3599, 7419, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 11119, 3, 6,
                                                                       10519, 7119, 10609, 3599,
                                                                       3659, 7519, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 11269, 3, 6,
                                                                       10699, 7239, 10789, 3779,
                                                                       3839, 7619, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 11419, 3, 6,
                                                                       10789, 7299, 10879, 3839,
                                                                       3899, 7719, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11569, 0, 6, 9919,
                                                                       6719, 9934, 4019, 4037,
                                                                       7819, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11614, 0, 6, 9934,
                                                                       6729, 9949, 4037, 4055,
                                                                       7849, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11659, 0, 6, 9949,
                                                                       6739, 9964, 4055, 4073,
                                                                       7879, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11704, 0, 6, 9964,
                                                                       6749, 9979, 4073, 4091,
                                                                       7909, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11749, 0, 6, 9994,
                                                                       6769, 10009, 4127, 4145,
                                                                       7939, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11794, 0, 6,
                                                                       10009, 6779, 10024, 4145,
                                                                       4163, 7969, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11839, 0, 6,
                                                                       10024, 6789, 10039, 4163,
                                                                       4181, 7999, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 11884, 0, 6,
                                                                       10039, 6799, 10054, 4181,
                                                                       4199, 8029, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 11929, 0, 3, 6,
                                                                       10069, 6819, 10114, 11569,
                                                                       7819, 11614, 4235, 4289,
                                                                       8059, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 12064, 0, 3, 6,
                                                                       10114, 6849, 10159, 11614,
                                                                       7849, 11659, 4289, 4343,
                                                                       8149, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 12199, 0, 3, 6,
                                                                       10159, 6879, 10204, 11659,
                                                                       7879, 11704, 4343, 4397,
                                                                       8239, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 12334, 0, 3, 6,
                                                                       10249, 6939, 10294, 11749,
                                                                       7939, 11794, 4505, 4559,
                                                                       8329, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 12469, 0, 3, 6,
                                                                       10294, 6969, 10339, 11794,
                                                                       7969, 11839, 4559, 4613,
                                                                       8419, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 12604, 0, 3, 6,
                                                                       10339, 6999, 10384, 11839,
                                                                       7999, 11884, 4613, 4667,
                                                                       8509, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 12739, 0, 3, 6,
                                                                       10429, 7059, 10519, 11929,
                                                                       8059, 12064, 4775, 4883,
                                                                       8599, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 13009, 0, 3, 6,
                                                                       10519, 7119, 10609, 12064,
                                                                       8149, 12199, 4883, 4991,
                                                                       8779, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 13279, 0, 3, 6,
                                                                       10699, 7239, 10789, 12334,
                                                                       8329, 12469, 5207, 5315,
                                                                       8959, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 13549, 0, 3, 6,
                                                                       10789, 7299, 10879, 12469,
                                                                       8419, 12604, 5315, 5423,
                                                                       9139, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 13819, 0, 3, 6,
                                                                       10969, 7419, 11119, 11929,
                                                                       12064, 12739, 8599, 13009,
                                                                       5639, 5819, 9319, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 14269, 0, 3, 6,
                                                                       11269, 7619, 11419, 12334,
                                                                       12469, 13279, 8959, 13549,
                                                                       6179, 6359, 9619, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 14719, 14269, 1, 150, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 14869, 14269, 1, 150, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 15019, 14269, 1, 150, ncols, alpha);

                    simdgeo::geom_s_x(buffer, 15169, 13819, 1, 150, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 15319, 13819, 1, 150, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 15469, 13819, 1, 150, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 15619, 14719, 900, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 16519, 15619, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 16519, 9, nmax);

        simdtrf::transform_g_inner(buffer, 16519, 15769, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 63 * nvalues + n * npairs, nvalues, buffer, 16519, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 16519, 15919, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 126 * nvalues + n * npairs, nvalues, buffer, 16519,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 16519, 16069, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 189 * nvalues + n * npairs, nvalues, buffer, 16519,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 16519, 16219, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 252 * nvalues + n * npairs, nvalues, buffer, 16519,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 16519, 16369, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 315 * nvalues + n * npairs, nvalues, buffer, 16519,
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
