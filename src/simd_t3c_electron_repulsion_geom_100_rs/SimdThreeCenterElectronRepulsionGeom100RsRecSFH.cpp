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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecSFH.hpp"

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
compute_rs_geom_100_sfh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_sfh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 281, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 284, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 287, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 290, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 293, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 296, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 299, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 302, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 305, 0, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 308, 0, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 311, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 314, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 317, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 320, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 323, 0, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 326, 0, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 329, 0, 6, 10, 11,
                                                                       29, 32, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 338, 0, 6, 11, 12,
                                                                       32, 35, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 347, 0, 6, 12, 13,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 356, 0, 6, 13, 14,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 365, 0, 6, 14, 15,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 374, 0, 6, 15, 16,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 383, 0, 6, 16, 17,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 392, 0, 6, 20, 21,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 401, 0, 6, 21, 22,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 410, 0, 6, 22, 23,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 419, 0, 6, 23, 24,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 428, 0, 6, 24, 25,
                                                                       65, 68, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 437, 0, 6, 25, 26,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 446, 0, 6, 26, 27,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 455, 0, 3, 6, 29,
                                                                       32, 77, 83, 329, 338,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 473, 0, 3, 6, 32,
                                                                       35, 83, 89, 338, 347,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 491, 0, 3, 6, 35,
                                                                       38, 89, 95, 347, 356,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 509, 0, 3, 6, 38,
                                                                       41, 95, 101, 356, 365,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 527, 0, 3, 6, 41,
                                                                       44, 101, 107, 365, 374,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 545, 0, 3, 6, 44,
                                                                       47, 107, 113, 374, 383,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 563, 0, 3, 6, 53,
                                                                       56, 119, 125, 392, 401,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 581, 0, 3, 6, 56,
                                                                       59, 125, 131, 401, 410,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 599, 0, 3, 6, 59,
                                                                       62, 131, 137, 410, 419,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 617, 0, 3, 6, 62,
                                                                       65, 137, 143, 419, 428,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 635, 0, 3, 6, 65,
                                                                       68, 143, 149, 428, 437,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 653, 0, 3, 6, 68,
                                                                       71, 149, 155, 437, 446,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 671, 0, 3, 6, 77,
                                                                       83, 161, 171, 455, 473,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 701, 0, 3, 6, 83,
                                                                       89, 171, 181, 473, 491,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 731, 0, 3, 6, 89,
                                                                       95, 181, 191, 491, 509,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 761, 0, 3, 6, 95,
                                                                       101, 191, 201, 509, 527,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 791, 0, 3, 6, 101,
                                                                       107, 201, 211, 527, 545,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 821, 0, 3, 6, 119,
                                                                       125, 221, 231, 563, 581,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 851, 0, 3, 6, 125,
                                                                       131, 231, 241, 581, 599,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 881, 0, 3, 6, 131,
                                                                       137, 241, 251, 599, 617,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 911, 0, 3, 6, 137,
                                                                       143, 251, 261, 617, 635,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 941, 0, 3, 6, 143,
                                                                       149, 261, 271, 635, 653,
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

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1133, 6, 29, 77,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1151, 6, 32, 83,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1169, 6, 35, 89,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1187, 6, 38, 95,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1205, 6, 41, 101,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1223, 6, 44, 107,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1241, 6, 47, 113,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1259, 6, 53, 119,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1277, 6, 56, 125,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1295, 6, 59, 131,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1313, 6, 62, 137,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1331, 6, 65, 143,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1349, 6, 68, 149,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1367, 6, 71, 155,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1385, 6, 77, 161,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1415, 6, 83, 171,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1445, 6, 89, 181,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1475, 6, 95, 191,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1505, 6, 101, 201,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1535, 6, 107, 211,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1565, 6, 119, 221,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1595, 6, 125, 231,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1625, 6, 131, 241,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1655, 6, 137, 251,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1685, 6, 143, 261,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1715, 6, 149, 271,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1745, 6, 10, 281,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1754, 6, 11, 284,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1763, 6, 12, 287,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1772, 6, 13, 290,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1781, 6, 14, 293,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1790, 6, 15, 296,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1799, 6, 16, 299,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1808, 6, 17, 302,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1817, 6, 20, 305,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1826, 6, 21, 308,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1835, 6, 22, 311,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1844, 6, 23, 314,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1853, 6, 24, 317,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1862, 6, 25, 320,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1871, 6, 26, 323,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1880, 6, 27, 326,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1889, 6, 29, 281,
                                                                       329, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1916, 6, 32, 284,
                                                                       338, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1943, 6, 35, 287,
                                                                       347, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1970, 6, 38, 290,
                                                                       356, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1997, 6, 41, 293,
                                                                       365, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2024, 6, 44, 296,
                                                                       374, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2051, 6, 47, 299,
                                                                       383, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2078, 6, 53, 305,
                                                                       392, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2105, 6, 56, 308,
                                                                       401, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2132, 6, 59, 311,
                                                                       410, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2159, 6, 62, 314,
                                                                       419, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2186, 6, 65, 317,
                                                                       428, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2213, 6, 68, 320,
                                                                       437, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2240, 6, 71, 323,
                                                                       446, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2267, 3, 6, 77,
                                                                       1889, 329, 1916, 455,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2321, 3, 6, 83,
                                                                       1916, 338, 1943, 473,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2375, 3, 6, 89,
                                                                       1943, 347, 1970, 491,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2429, 3, 6, 95,
                                                                       1970, 356, 1997, 509,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2483, 3, 6, 101,
                                                                       1997, 365, 2024, 527,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2537, 3, 6, 107,
                                                                       2024, 374, 2051, 545,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2591, 3, 6, 119,
                                                                       2078, 392, 2105, 563,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2645, 3, 6, 125,
                                                                       2105, 401, 2132, 581,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2699, 3, 6, 131,
                                                                       2132, 410, 2159, 599,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2753, 3, 6, 137,
                                                                       2159, 419, 2186, 617,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2807, 3, 6, 143,
                                                                       2186, 428, 2213, 635,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2861, 3, 6, 149,
                                                                       2213, 437, 2240, 653,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 2915, 3, 6, 161,
                                                                       2267, 455, 2321, 671,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3005, 3, 6, 171,
                                                                       2321, 473, 2375, 701,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3095, 3, 6, 181,
                                                                       2375, 491, 2429, 731,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3185, 3, 6, 191,
                                                                       2429, 509, 2483, 761,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3275, 3, 6, 201,
                                                                       2483, 527, 2537, 791,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3365, 3, 6, 221,
                                                                       2591, 563, 2645, 821,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3455, 3, 6, 231,
                                                                       2645, 581, 2699, 851,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3545, 3, 6, 241,
                                                                       2699, 599, 2753, 881,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3635, 3, 6, 251,
                                                                       2753, 617, 2807, 911,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3725, 3, 6, 261,
                                                                       2807, 635, 2861, 941,
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

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4115, 3, 6, 3899,
                                                                       1025, 3917, 77, 83, 1169,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4151, 3, 6, 3917,
                                                                       1034, 3935, 83, 89, 1187,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4187, 3, 6, 3935,
                                                                       1043, 3953, 89, 95, 1205,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4223, 3, 6, 3953,
                                                                       1052, 3971, 95, 101, 1223,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4259, 3, 6, 3971,
                                                                       1061, 3989, 101, 107,
                                                                       1241, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4295, 3, 6, 4007,
                                                                       1079, 4025, 119, 125,
                                                                       1295, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4331, 3, 6, 4025,
                                                                       1088, 4043, 125, 131,
                                                                       1313, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4367, 3, 6, 4043,
                                                                       1097, 4061, 131, 137,
                                                                       1331, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4403, 3, 6, 4061,
                                                                       1106, 4079, 137, 143,
                                                                       1349, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4439, 3, 6, 4079,
                                                                       1115, 4097, 143, 149,
                                                                       1367, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4475, 3, 6, 4115,
                                                                       1169, 4151, 161, 171,
                                                                       1445, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4535, 3, 6, 4151,
                                                                       1187, 4187, 171, 181,
                                                                       1475, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4595, 3, 6, 4187,
                                                                       1205, 4223, 181, 191,
                                                                       1505, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4655, 3, 6, 4223,
                                                                       1223, 4259, 191, 201,
                                                                       1535, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4715, 3, 6, 4295,
                                                                       1295, 4331, 221, 231,
                                                                       1625, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4775, 3, 6, 4331,
                                                                       1313, 4367, 231, 241,
                                                                       1655, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4835, 3, 6, 4367,
                                                                       1331, 4403, 241, 251,
                                                                       1685, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 4895, 3, 6, 4403,
                                                                       1349, 4439, 251, 261,
                                                                       1715, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4955, 0, 6, 3815,
                                                                       977, 3821, 1763, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4973, 0, 6, 3821,
                                                                       980, 3827, 1772, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4991, 0, 6, 3827,
                                                                       983, 3833, 1781, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5009, 0, 6, 3833,
                                                                       986, 3839, 1790, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5027, 0, 6, 3839,
                                                                       989, 3845, 1799, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5045, 0, 6, 3845,
                                                                       992, 3851, 1808, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5063, 0, 6, 3857,
                                                                       1004, 3863, 1835, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5081, 0, 6, 3863,
                                                                       1007, 3869, 1844, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5099, 0, 6, 3869,
                                                                       1010, 3875, 1853, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5117, 0, 6, 3875,
                                                                       1013, 3881, 1862, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5135, 0, 6, 3881,
                                                                       1016, 3887, 1871, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5153, 0, 6, 3887,
                                                                       1019, 3893, 1880, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5171, 0, 3, 6,
                                                                       3899, 1025, 3917, 4955,
                                                                       1763, 4973, 329, 338,
                                                                       1943, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5225, 0, 3, 6,
                                                                       3917, 1034, 3935, 4973,
                                                                       1772, 4991, 338, 347,
                                                                       1970, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5279, 0, 3, 6,
                                                                       3935, 1043, 3953, 4991,
                                                                       1781, 5009, 347, 356,
                                                                       1997, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5333, 0, 3, 6,
                                                                       3953, 1052, 3971, 5009,
                                                                       1790, 5027, 356, 365,
                                                                       2024, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5387, 0, 3, 6,
                                                                       3971, 1061, 3989, 5027,
                                                                       1799, 5045, 365, 374,
                                                                       2051, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5441, 0, 3, 6,
                                                                       4007, 1079, 4025, 5063,
                                                                       1835, 5081, 392, 401,
                                                                       2132, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5495, 0, 3, 6,
                                                                       4025, 1088, 4043, 5081,
                                                                       1844, 5099, 401, 410,
                                                                       2159, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5549, 0, 3, 6,
                                                                       4043, 1097, 4061, 5099,
                                                                       1853, 5117, 410, 419,
                                                                       2186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5603, 0, 3, 6,
                                                                       4061, 1106, 4079, 5117,
                                                                       1862, 5135, 419, 428,
                                                                       2213, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5657, 0, 3, 6,
                                                                       4079, 1115, 4097, 5135,
                                                                       1871, 5153, 428, 437,
                                                                       2240, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5711, 0, 3, 6,
                                                                       4115, 1169, 4151, 5171,
                                                                       1943, 5225, 455, 473,
                                                                       2375, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5819, 0, 3, 6,
                                                                       4151, 1187, 4187, 5225,
                                                                       1970, 5279, 473, 491,
                                                                       2429, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 5927, 0, 3, 6,
                                                                       4187, 1205, 4223, 5279,
                                                                       1997, 5333, 491, 509,
                                                                       2483, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6035, 0, 3, 6,
                                                                       4223, 1223, 4259, 5333,
                                                                       2024, 5387, 509, 527,
                                                                       2537, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6143, 0, 3, 6,
                                                                       4295, 1295, 4331, 5441,
                                                                       2132, 5495, 563, 581,
                                                                       2699, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6251, 0, 3, 6,
                                                                       4331, 1313, 4367, 5495,
                                                                       2159, 5549, 581, 599,
                                                                       2753, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6359, 0, 3, 6,
                                                                       4367, 1331, 4403, 5549,
                                                                       2186, 5603, 599, 617,
                                                                       2807, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6467, 0, 3, 6,
                                                                       4403, 1349, 4439, 5603,
                                                                       2213, 5657, 617, 635,
                                                                       2861, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 6575, 0, 3, 6,
                                                                       4475, 1445, 4535, 5171,
                                                                       5225, 5711, 2375, 5819,
                                                                       671, 701, 3095, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 6755, 0, 3, 6,
                                                                       4535, 1475, 4595, 5225,
                                                                       5279, 5819, 2429, 5927,
                                                                       701, 731, 3185, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 6935, 0, 3, 6,
                                                                       4595, 1505, 4655, 5279,
                                                                       5333, 5927, 2483, 6035,
                                                                       731, 761, 3275, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 7115, 0, 3, 6,
                                                                       4715, 1625, 4775, 5441,
                                                                       5495, 6143, 2699, 6251,
                                                                       821, 851, 3545, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 7295, 0, 3, 6,
                                                                       4775, 1655, 4835, 5495,
                                                                       5549, 6251, 2753, 6359,
                                                                       851, 881, 3635, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 7475, 0, 3, 6,
                                                                       4835, 1685, 4895, 5549,
                                                                       5603, 6359, 2807, 6467,
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

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8155, 3, 6, 7795,
                                                                       3899, 7825, 1133, 1151,
                                                                       4115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8215, 3, 6, 7825,
                                                                       3917, 7855, 1151, 1169,
                                                                       4151, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8275, 3, 6, 7855,
                                                                       3935, 7885, 1169, 1187,
                                                                       4187, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8335, 3, 6, 7885,
                                                                       3953, 7915, 1187, 1205,
                                                                       4223, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8395, 3, 6, 7915,
                                                                       3971, 7945, 1205, 1223,
                                                                       4259, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8455, 3, 6, 7975,
                                                                       4007, 8005, 1259, 1277,
                                                                       4295, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8515, 3, 6, 8005,
                                                                       4025, 8035, 1277, 1295,
                                                                       4331, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8575, 3, 6, 8035,
                                                                       4043, 8065, 1295, 1313,
                                                                       4367, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8635, 3, 6, 8065,
                                                                       4061, 8095, 1313, 1331,
                                                                       4403, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8695, 3, 6, 8095,
                                                                       4079, 8125, 1331, 1349,
                                                                       4439, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 8755, 3, 6, 8155,
                                                                       4115, 8215, 1385, 1415,
                                                                       4475, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 8855, 3, 6, 8215,
                                                                       4151, 8275, 1415, 1445,
                                                                       4535, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 8955, 3, 6, 8275,
                                                                       4187, 8335, 1445, 1475,
                                                                       4595, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9055, 3, 6, 8335,
                                                                       4223, 8395, 1475, 1505,
                                                                       4655, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9155, 3, 6, 8455,
                                                                       4295, 8515, 1565, 1595,
                                                                       4715, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9255, 3, 6, 8515,
                                                                       4331, 8575, 1595, 1625,
                                                                       4775, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9355, 3, 6, 8575,
                                                                       4367, 8635, 1625, 1655,
                                                                       4835, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9455, 3, 6, 8635,
                                                                       4403, 8695, 1655, 1685,
                                                                       4895, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9555, 0, 6, 7655,
                                                                       3815, 7665, 1745, 1754,
                                                                       4955, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9585, 0, 6, 7665,
                                                                       3821, 7675, 1754, 1763,
                                                                       4973, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9615, 0, 6, 7675,
                                                                       3827, 7685, 1763, 1772,
                                                                       4991, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9645, 0, 6, 7685,
                                                                       3833, 7695, 1772, 1781,
                                                                       5009, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9675, 0, 6, 7695,
                                                                       3839, 7705, 1781, 1790,
                                                                       5027, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9705, 0, 6, 7705,
                                                                       3845, 7715, 1790, 1799,
                                                                       5045, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9735, 0, 6, 7725,
                                                                       3857, 7735, 1817, 1826,
                                                                       5063, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9765, 0, 6, 7735,
                                                                       3863, 7745, 1826, 1835,
                                                                       5081, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9795, 0, 6, 7745,
                                                                       3869, 7755, 1835, 1844,
                                                                       5099, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9825, 0, 6, 7755,
                                                                       3875, 7765, 1844, 1853,
                                                                       5117, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9855, 0, 6, 7765,
                                                                       3881, 7775, 1853, 1862,
                                                                       5135, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9885, 0, 6, 7775,
                                                                       3887, 7785, 1862, 1871,
                                                                       5153, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9915, 0, 3, 6,
                                                                       7795, 3899, 7825, 9555,
                                                                       4955, 9585, 1889, 1916,
                                                                       5171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10005, 0, 3, 6,
                                                                       7825, 3917, 7855, 9585,
                                                                       4973, 9615, 1916, 1943,
                                                                       5225, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10095, 0, 3, 6,
                                                                       7855, 3935, 7885, 9615,
                                                                       4991, 9645, 1943, 1970,
                                                                       5279, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10185, 0, 3, 6,
                                                                       7885, 3953, 7915, 9645,
                                                                       5009, 9675, 1970, 1997,
                                                                       5333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10275, 0, 3, 6,
                                                                       7915, 3971, 7945, 9675,
                                                                       5027, 9705, 1997, 2024,
                                                                       5387, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10365, 0, 3, 6,
                                                                       7975, 4007, 8005, 9735,
                                                                       5063, 9765, 2078, 2105,
                                                                       5441, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10455, 0, 3, 6,
                                                                       8005, 4025, 8035, 9765,
                                                                       5081, 9795, 2105, 2132,
                                                                       5495, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10545, 0, 3, 6,
                                                                       8035, 4043, 8065, 9795,
                                                                       5099, 9825, 2132, 2159,
                                                                       5549, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10635, 0, 3, 6,
                                                                       8065, 4061, 8095, 9825,
                                                                       5117, 9855, 2159, 2186,
                                                                       5603, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10725, 0, 3, 6,
                                                                       8095, 4079, 8125, 9855,
                                                                       5135, 9885, 2186, 2213,
                                                                       5657, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 10815, 0, 3, 6,
                                                                       8155, 4115, 8215, 9915,
                                                                       5171, 10005, 2267, 2321,
                                                                       5711, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 10995, 0, 3, 6,
                                                                       8215, 4151, 8275, 10005,
                                                                       5225, 10095, 2321, 2375,
                                                                       5819, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 11175, 0, 3, 6,
                                                                       8275, 4187, 8335, 10095,
                                                                       5279, 10185, 2375, 2429,
                                                                       5927, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 11355, 0, 3, 6,
                                                                       8335, 4223, 8395, 10185,
                                                                       5333, 10275, 2429, 2483,
                                                                       6035, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 11535, 0, 3, 6,
                                                                       8455, 4295, 8515, 10365,
                                                                       5441, 10455, 2591, 2645,
                                                                       6143, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 11715, 0, 3, 6,
                                                                       8515, 4331, 8575, 10455,
                                                                       5495, 10545, 2645, 2699,
                                                                       6251, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 11895, 0, 3, 6,
                                                                       8575, 4367, 8635, 10545,
                                                                       5549, 10635, 2699, 2753,
                                                                       6359, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 12075, 0, 3, 6,
                                                                       8635, 4403, 8695, 10635,
                                                                       5603, 10725, 2753, 2807,
                                                                       6467, ncols, gamma, p,
                                                                       q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 12255, 0, 3, 6,
                                                                       8755, 4475, 8855, 9915,
                                                                       10005, 10815, 5711, 10995,
                                                                       2915, 3005, 6575, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 12555, 0, 3, 6,
                                                                       8855, 4535, 8955, 10005,
                                                                       10095, 10995, 5819, 11175,
                                                                       3005, 3095, 6755, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 12855, 0, 3, 6,
                                                                       8955, 4595, 9055, 10095,
                                                                       10185, 11175, 5927, 11355,
                                                                       3095, 3185, 6935, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 13155, 0, 3, 6,
                                                                       9155, 4715, 9255, 10365,
                                                                       10455, 11535, 6143, 11715,
                                                                       3365, 3455, 7115, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 13455, 0, 3, 6,
                                                                       9255, 4775, 9355, 10455,
                                                                       10545, 11715, 6251, 11895,
                                                                       3455, 3545, 7295, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 13755, 0, 3, 6,
                                                                       9355, 4835, 9455, 10545,
                                                                       10635, 11895, 6359, 12075,
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

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14565, 3, 6,
                                                                       14205, 7855, 14250, 4115,
                                                                       4151, 8275, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14655, 3, 6,
                                                                       14250, 7885, 14295, 4151,
                                                                       4187, 8335, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14745, 3, 6,
                                                                       14295, 7915, 14340, 4187,
                                                                       4223, 8395, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14835, 3, 6,
                                                                       14385, 8035, 14430, 4295,
                                                                       4331, 8575, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 14925, 3, 6,
                                                                       14430, 8065, 14475, 4331,
                                                                       4367, 8635, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 15015, 3, 6,
                                                                       14475, 8095, 14520, 4367,
                                                                       4403, 8695, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 15105, 3, 6,
                                                                       14565, 8275, 14655, 4475,
                                                                       4535, 8955, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 15255, 3, 6,
                                                                       14655, 8335, 14745, 4535,
                                                                       4595, 9055, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 15405, 3, 6,
                                                                       14835, 8575, 14925, 4715,
                                                                       4775, 9355, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 15555, 3, 6,
                                                                       14925, 8635, 15015, 4775,
                                                                       4835, 9455, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15705, 0, 6,
                                                                       14055, 7675, 14070, 4955,
                                                                       4973, 9615, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15750, 0, 6,
                                                                       14070, 7685, 14085, 4973,
                                                                       4991, 9645, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15795, 0, 6,
                                                                       14085, 7695, 14100, 4991,
                                                                       5009, 9675, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15840, 0, 6,
                                                                       14100, 7705, 14115, 5009,
                                                                       5027, 9705, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15885, 0, 6,
                                                                       14130, 7745, 14145, 5063,
                                                                       5081, 9795, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15930, 0, 6,
                                                                       14145, 7755, 14160, 5081,
                                                                       5099, 9825, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15975, 0, 6,
                                                                       14160, 7765, 14175, 5099,
                                                                       5117, 9855, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 16020, 0, 6,
                                                                       14175, 7775, 14190, 5117,
                                                                       5135, 9885, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 16065, 0, 3, 6,
                                                                       14205, 7855, 14250, 15705,
                                                                       9615, 15750, 5171, 5225,
                                                                       10095, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 16200, 0, 3, 6,
                                                                       14250, 7885, 14295, 15750,
                                                                       9645, 15795, 5225, 5279,
                                                                       10185, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 16335, 0, 3, 6,
                                                                       14295, 7915, 14340, 15795,
                                                                       9675, 15840, 5279, 5333,
                                                                       10275, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 16470, 0, 3, 6,
                                                                       14385, 8035, 14430, 15885,
                                                                       9795, 15930, 5441, 5495,
                                                                       10545, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 16605, 0, 3, 6,
                                                                       14430, 8065, 14475, 15930,
                                                                       9825, 15975, 5495, 5549,
                                                                       10635, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 16740, 0, 3, 6,
                                                                       14475, 8095, 14520, 15975,
                                                                       9855, 16020, 5549, 5603,
                                                                       10725, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 16875, 0, 3, 6,
                                                                       14565, 8275, 14655, 16065,
                                                                       10095, 16200, 5711, 5819,
                                                                       11175, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 17145, 0, 3, 6,
                                                                       14655, 8335, 14745, 16200,
                                                                       10185, 16335, 5819, 5927,
                                                                       11355, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 17415, 0, 3, 6,
                                                                       14835, 8575, 14925, 16470,
                                                                       10545, 16605, 6143, 6251,
                                                                       11895, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 17685, 0, 3, 6,
                                                                       14925, 8635, 15015, 16605,
                                                                       10635, 16740, 6251, 6359,
                                                                       12075, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 17955, 0, 3, 6,
                                                                       15105, 8955, 15255, 16065,
                                                                       16200, 16875, 11175,
                                                                       17145, 6575, 6755, 12855,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 18405, 0, 3, 6,
                                                                       15405, 9355, 15555, 16470,
                                                                       16605, 17415, 11895,
                                                                       17685, 7115, 7295, 13755,
                                                                       ncols, gamma, p, q);

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

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 19569, 3, 6,
                                                                       19065, 14205, 19128, 8155,
                                                                       8215, 14565, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 19695, 3, 6,
                                                                       19128, 14250, 19191, 8215,
                                                                       8275, 14655, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 19821, 3, 6,
                                                                       19191, 14295, 19254, 8275,
                                                                       8335, 14745, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 19947, 3, 6,
                                                                       19317, 14385, 19380, 8455,
                                                                       8515, 14835, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 20073, 3, 6,
                                                                       19380, 14430, 19443, 8515,
                                                                       8575, 14925, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 20199, 3, 6,
                                                                       19443, 14475, 19506, 8575,
                                                                       8635, 15015, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 20325, 3, 6,
                                                                       19569, 14565, 19695, 8755,
                                                                       8855, 15105, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 20535, 3, 6,
                                                                       19695, 14655, 19821, 8855,
                                                                       8955, 15255, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 20745, 3, 6,
                                                                       19947, 14835, 20073, 9155,
                                                                       9255, 15405, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 20955, 3, 6,
                                                                       20073, 14925, 20199, 9255,
                                                                       9355, 15555, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21165, 0, 6,
                                                                       18855, 14055, 18876, 9555,
                                                                       9585, 15705, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21228, 0, 6,
                                                                       18876, 14070, 18897, 9585,
                                                                       9615, 15750, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21291, 0, 6,
                                                                       18897, 14085, 18918, 9615,
                                                                       9645, 15795, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21354, 0, 6,
                                                                       18918, 14100, 18939, 9645,
                                                                       9675, 15840, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21417, 0, 6,
                                                                       18960, 14130, 18981, 9735,
                                                                       9765, 15885, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21480, 0, 6,
                                                                       18981, 14145, 19002, 9765,
                                                                       9795, 15930, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21543, 0, 6,
                                                                       19002, 14160, 19023, 9795,
                                                                       9825, 15975, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 21606, 0, 6,
                                                                       19023, 14175, 19044, 9825,
                                                                       9855, 16020, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 21669, 0, 3, 6,
                                                                       19065, 14205, 19128,
                                                                       21165, 15705, 21228, 9915,
                                                                       10005, 16065, ncols,
                                                                       gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 21858, 0, 3, 6,
                                                                       19128, 14250, 19191,
                                                                       21228, 15750, 21291,
                                                                       10005, 10095, 16200,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 22047, 0, 3, 6,
                                                                       19191, 14295, 19254,
                                                                       21291, 15795, 21354,
                                                                       10095, 10185, 16335,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 22236, 0, 3, 6,
                                                                       19317, 14385, 19380,
                                                                       21417, 15885, 21480,
                                                                       10365, 10455, 16470,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 22425, 0, 3, 6,
                                                                       19380, 14430, 19443,
                                                                       21480, 15930, 21543,
                                                                       10455, 10545, 16605,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 22614, 0, 3, 6,
                                                                       19443, 14475, 19506,
                                                                       21543, 15975, 21606,
                                                                       10545, 10635, 16740,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 22803, 0, 3, 6,
                                                                       19569, 14565, 19695,
                                                                       21669, 16065, 21858,
                                                                       10815, 10995, 16875,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 23181, 0, 3, 6,
                                                                       19695, 14655, 19821,
                                                                       21858, 16200, 22047,
                                                                       10995, 11175, 17145,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 23559, 0, 3, 6,
                                                                       19947, 14835, 20073,
                                                                       22236, 16470, 22425,
                                                                       11535, 11715, 17415,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 23937, 0, 3, 6,
                                                                       20073, 14925, 20199,
                                                                       22425, 16605, 22614,
                                                                       11715, 11895, 17685,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 24315, 0, 3, 6,
                                                                       20325, 15105, 20535,
                                                                       21669, 21858, 22803,
                                                                       16875, 23181, 12255,
                                                                       12555, 17955, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 24945, 0, 3, 6,
                                                                       20745, 15405, 20955,
                                                                       22236, 22425, 23559,
                                                                       17415, 23937, 13155,
                                                                       13455, 18405, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 25575, 24945, 1, 210, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 25785, 24945, 1, 210, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 25995, 24945, 1, 210, ncols, alpha);

                    simdgeo::geom_s_x(buffer, 26205, 24315, 1, 210, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 26415, 24315, 1, 210, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 26625, 24315, 1, 210, ncols, alpha);

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
