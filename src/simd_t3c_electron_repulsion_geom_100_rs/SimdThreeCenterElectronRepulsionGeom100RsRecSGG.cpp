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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecSGG.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGS.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
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
compute_rs_geom_100_sgg_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_sgg_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 31044, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 31044, 29559, 1350, dimensions);

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

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 85, 3, 6, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 91, 3, 6, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 97, 3, 6, 12, 13,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 103, 3, 6, 13, 14,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 109, 3, 6, 14, 15,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 115, 3, 6, 15, 16,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 121, 3, 6, 16, 17,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 127, 3, 6, 17, 18,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 133, 3, 6, 21, 22,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 139, 3, 6, 22, 23,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 145, 3, 6, 23, 24,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 151, 3, 6, 24, 25,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 157, 3, 6, 25, 26,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 163, 3, 6, 26, 27,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 169, 3, 6, 27, 28,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 175, 3, 6, 28, 29,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 181, 3, 6, 31, 34,
                                                                       85, 91, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 191, 3, 6, 34, 37,
                                                                       91, 97, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 201, 3, 6, 37, 40,
                                                                       97, 103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 211, 3, 6, 40, 43,
                                                                       103, 109, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 221, 3, 6, 43, 46,
                                                                       109, 115, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 231, 3, 6, 46, 49,
                                                                       115, 121, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 241, 3, 6, 49, 52,
                                                                       121, 127, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 251, 3, 6, 58, 61,
                                                                       133, 139, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 261, 3, 6, 61, 64,
                                                                       139, 145, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 271, 3, 6, 64, 67,
                                                                       145, 151, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 281, 3, 6, 67, 70,
                                                                       151, 157, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 291, 3, 6, 70, 73,
                                                                       157, 163, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 301, 3, 6, 73, 76,
                                                                       163, 169, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 311, 3, 6, 76, 79,
                                                                       169, 175, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 321, 3, 6, 85, 91,
                                                                       181, 191, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 336, 3, 6, 91, 97,
                                                                       191, 201, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 351, 3, 6, 97,
                                                                       103, 201, 211, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 366, 3, 6, 103,
                                                                       109, 211, 221, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 381, 3, 6, 109,
                                                                       115, 221, 231, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 396, 3, 6, 115,
                                                                       121, 231, 241, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 411, 3, 6, 133,
                                                                       139, 251, 261, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 426, 3, 6, 139,
                                                                       145, 261, 271, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 441, 3, 6, 145,
                                                                       151, 271, 281, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 456, 3, 6, 151,
                                                                       157, 281, 291, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 471, 3, 6, 157,
                                                                       163, 291, 301, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 486, 3, 6, 163,
                                                                       169, 301, 311, ncols,
                                                                       gamma, p, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 501, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 504, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 507, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 510, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 513, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 516, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 519, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 522, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 525, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 528, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 531, 0, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 534, 0, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 537, 0, 6, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 540, 0, 6, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 543, 0, 6, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 552, 0, 6, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 561, 0, 6, 12, 13,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 570, 0, 6, 13, 14,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 579, 0, 6, 14, 15,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 588, 0, 6, 15, 16,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 597, 0, 6, 16, 17,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 606, 0, 6, 17, 18,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 615, 0, 6, 21, 22,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 624, 0, 6, 22, 23,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 633, 0, 6, 23, 24,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 642, 0, 6, 24, 25,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 651, 0, 6, 25, 26,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 660, 0, 6, 26, 27,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 669, 0, 6, 27, 28,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 678, 0, 6, 28, 29,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 687, 0, 3, 6, 31,
                                                                       34, 85, 91, 543, 552,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 705, 0, 3, 6, 34,
                                                                       37, 91, 97, 552, 561,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 723, 0, 3, 6, 37,
                                                                       40, 97, 103, 561, 570,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 741, 0, 3, 6, 40,
                                                                       43, 103, 109, 570, 579,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 759, 0, 3, 6, 43,
                                                                       46, 109, 115, 579, 588,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 777, 0, 3, 6, 46,
                                                                       49, 115, 121, 588, 597,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 795, 0, 3, 6, 49,
                                                                       52, 121, 127, 597, 606,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 813, 0, 3, 6, 58,
                                                                       61, 133, 139, 615, 624,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 831, 0, 3, 6, 61,
                                                                       64, 139, 145, 624, 633,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 849, 0, 3, 6, 64,
                                                                       67, 145, 151, 633, 642,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 867, 0, 3, 6, 67,
                                                                       70, 151, 157, 642, 651,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 885, 0, 3, 6, 70,
                                                                       73, 157, 163, 651, 660,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 903, 0, 3, 6, 73,
                                                                       76, 163, 169, 660, 669,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 921, 0, 3, 6, 76,
                                                                       79, 169, 175, 669, 678,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 939, 0, 3, 6, 85,
                                                                       91, 181, 191, 687, 705,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 969, 0, 3, 6, 91,
                                                                       97, 191, 201, 705, 723,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 999, 0, 3, 6, 97,
                                                                       103, 201, 211, 723, 741,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1029, 0, 3, 6,
                                                                       103, 109, 211, 221, 741,
                                                                       759, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1059, 0, 3, 6,
                                                                       109, 115, 221, 231, 759,
                                                                       777, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1089, 0, 3, 6,
                                                                       115, 121, 231, 241, 777,
                                                                       795, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1119, 0, 3, 6,
                                                                       133, 139, 251, 261, 813,
                                                                       831, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1149, 0, 3, 6,
                                                                       139, 145, 261, 271, 831,
                                                                       849, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1179, 0, 3, 6,
                                                                       145, 151, 271, 281, 849,
                                                                       867, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1209, 0, 3, 6,
                                                                       151, 157, 281, 291, 867,
                                                                       885, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1239, 0, 3, 6,
                                                                       157, 163, 291, 301, 885,
                                                                       903, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1269, 0, 3, 6,
                                                                       163, 169, 301, 311, 903,
                                                                       921, ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1299, 0, 3, 6,
                                                                       181, 191, 321, 336, 939,
                                                                       969, ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1344, 0, 3, 6,
                                                                       191, 201, 336, 351, 969,
                                                                       999, ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1389, 0, 3, 6,
                                                                       201, 211, 351, 366, 999,
                                                                       1029, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1434, 0, 3, 6,
                                                                       211, 221, 366, 381, 1029,
                                                                       1059, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1479, 0, 3, 6,
                                                                       221, 231, 381, 396, 1059,
                                                                       1089, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1524, 0, 3, 6,
                                                                       251, 261, 411, 426, 1119,
                                                                       1149, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1569, 0, 3, 6,
                                                                       261, 271, 426, 441, 1149,
                                                                       1179, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1614, 0, 3, 6,
                                                                       271, 281, 441, 456, 1179,
                                                                       1209, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1659, 0, 3, 6,
                                                                       281, 291, 456, 471, 1209,
                                                                       1239, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1704, 0, 3, 6,
                                                                       291, 301, 471, 486, 1239,
                                                                       1269, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1749, 6, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1752, 6, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1755, 6, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1758, 6, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1761, 6, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1764, 6, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1767, 6, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1770, 6, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1773, 6, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1776, 6, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1779, 6, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1782, 6, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1785, 6, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1788, 6, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1791, 6, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1794, 6, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1797, 6, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1806, 6, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1815, 6, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1824, 6, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1833, 6, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1842, 6, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1851, 6, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1860, 6, 23, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1869, 6, 24, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1878, 6, 25, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1887, 6, 26, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1896, 6, 27, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1905, 6, 28, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1914, 6, 29, 82,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1923, 6, 37, 97,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1941, 6, 40, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1959, 6, 43, 109,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1977, 6, 46, 115,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1995, 6, 49, 121,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2013, 6, 52, 127,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2031, 6, 64, 145,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2049, 6, 67, 151,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2067, 6, 70, 157,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2085, 6, 73, 163,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2103, 6, 76, 169,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2121, 6, 79, 175,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2139, 6, 97, 201,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2169, 6, 103, 211,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2199, 6, 109, 221,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2229, 6, 115, 231,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2259, 6, 121, 241,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2289, 6, 145, 271,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2319, 6, 151, 281,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2349, 6, 157, 291,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2379, 6, 163, 301,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2409, 6, 169, 311,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2439, 6, 201, 351,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2484, 6, 211, 366,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2529, 6, 221, 381,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2574, 6, 231, 396,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2619, 6, 271, 441,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2664, 6, 281, 456,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2709, 6, 291, 471,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2754, 6, 301, 486,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2799, 6, 12, 501,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2808, 6, 13, 504,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2817, 6, 14, 507,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2826, 6, 15, 510,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2835, 6, 16, 513,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2844, 6, 17, 516,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2853, 6, 18, 519,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2862, 6, 23, 522,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2871, 6, 24, 525,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2880, 6, 25, 528,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2889, 6, 26, 531,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2898, 6, 27, 534,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2907, 6, 28, 537,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2916, 6, 29, 540,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2925, 6, 37, 501,
                                                                       561, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2952, 6, 40, 504,
                                                                       570, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2979, 6, 43, 507,
                                                                       579, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3006, 6, 46, 510,
                                                                       588, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3033, 6, 49, 513,
                                                                       597, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3060, 6, 52, 516,
                                                                       606, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3087, 6, 64, 522,
                                                                       633, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3114, 6, 67, 525,
                                                                       642, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3141, 6, 70, 528,
                                                                       651, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3168, 6, 73, 531,
                                                                       660, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3195, 6, 76, 534,
                                                                       669, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3222, 6, 79, 537,
                                                                       678, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3249, 3, 6, 97,
                                                                       2925, 561, 2952, 723,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3303, 3, 6, 103,
                                                                       2952, 570, 2979, 741,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3357, 3, 6, 109,
                                                                       2979, 579, 3006, 759,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3411, 3, 6, 115,
                                                                       3006, 588, 3033, 777,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3465, 3, 6, 121,
                                                                       3033, 597, 3060, 795,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3519, 3, 6, 145,
                                                                       3087, 633, 3114, 849,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3573, 3, 6, 151,
                                                                       3114, 642, 3141, 867,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3627, 3, 6, 157,
                                                                       3141, 651, 3168, 885,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3681, 3, 6, 163,
                                                                       3168, 660, 3195, 903,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3735, 3, 6, 169,
                                                                       3195, 669, 3222, 921,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3789, 3, 6, 201,
                                                                       3249, 723, 3303, 999,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3879, 3, 6, 211,
                                                                       3303, 741, 3357, 1029,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3969, 3, 6, 221,
                                                                       3357, 759, 3411, 1059,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 4059, 3, 6, 231,
                                                                       3411, 777, 3465, 1089,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 4149, 3, 6, 271,
                                                                       3519, 849, 3573, 1179,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 4239, 3, 6, 281,
                                                                       3573, 867, 3627, 1209,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 4329, 3, 6, 291,
                                                                       3627, 885, 3681, 1239,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 4419, 3, 6, 301,
                                                                       3681, 903, 3735, 1269,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 4509, 3, 6, 351,
                                                                       3789, 999, 3879, 1389,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 4644, 3, 6, 366,
                                                                       3879, 1029, 3969, 1434,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 4779, 3, 6, 381,
                                                                       3969, 1059, 4059, 1479,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 4914, 3, 6, 441,
                                                                       4149, 1179, 4239, 1614,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 5049, 3, 6, 456,
                                                                       4239, 1209, 4329, 1659,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 5184, 3, 6, 471,
                                                                       4329, 1239, 4419, 1704,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5319, 6, 10, 11,
                                                                       1749, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5325, 6, 11, 12,
                                                                       1752, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5331, 6, 12, 13,
                                                                       1755, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5337, 6, 13, 14,
                                                                       1758, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5343, 6, 14, 15,
                                                                       1761, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5349, 6, 15, 16,
                                                                       1764, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5355, 6, 16, 17,
                                                                       1767, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5361, 6, 17, 18,
                                                                       1770, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5367, 6, 21, 22,
                                                                       1773, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5373, 6, 22, 23,
                                                                       1776, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5379, 6, 23, 24,
                                                                       1779, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5385, 6, 24, 25,
                                                                       1782, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5391, 6, 25, 26,
                                                                       1785, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5397, 6, 26, 27,
                                                                       1788, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5403, 6, 27, 28,
                                                                       1791, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 5409, 6, 28, 29,
                                                                       1794, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5415, 3, 6, 5319,
                                                                       1749, 5325, 1797, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5433, 3, 6, 5325,
                                                                       1752, 5331, 1806, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5451, 3, 6, 5331,
                                                                       1755, 5337, 1815, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5469, 3, 6, 5337,
                                                                       1758, 5343, 1824, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5487, 3, 6, 5343,
                                                                       1761, 5349, 1833, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5505, 3, 6, 5349,
                                                                       1764, 5355, 1842, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5523, 3, 6, 5355,
                                                                       1767, 5361, 1851, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5541, 3, 6, 5367,
                                                                       1773, 5373, 1860, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5559, 3, 6, 5373,
                                                                       1776, 5379, 1869, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5577, 3, 6, 5379,
                                                                       1779, 5385, 1878, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5595, 3, 6, 5385,
                                                                       1782, 5391, 1887, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5613, 3, 6, 5391,
                                                                       1785, 5397, 1896, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5631, 3, 6, 5397,
                                                                       1788, 5403, 1905, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 5649, 3, 6, 5403,
                                                                       1791, 5409, 1914, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5667, 3, 6, 5415,
                                                                       1797, 5433, 85, 91, 1923,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5703, 3, 6, 5433,
                                                                       1806, 5451, 91, 97, 1941,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5739, 3, 6, 5451,
                                                                       1815, 5469, 97, 103, 1959,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5775, 3, 6, 5469,
                                                                       1824, 5487, 103, 109,
                                                                       1977, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5811, 3, 6, 5487,
                                                                       1833, 5505, 109, 115,
                                                                       1995, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5847, 3, 6, 5505,
                                                                       1842, 5523, 115, 121,
                                                                       2013, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5883, 3, 6, 5541,
                                                                       1860, 5559, 133, 139,
                                                                       2031, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5919, 3, 6, 5559,
                                                                       1869, 5577, 139, 145,
                                                                       2049, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5955, 3, 6, 5577,
                                                                       1878, 5595, 145, 151,
                                                                       2067, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5991, 3, 6, 5595,
                                                                       1887, 5613, 151, 157,
                                                                       2085, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6027, 3, 6, 5613,
                                                                       1896, 5631, 157, 163,
                                                                       2103, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 6063, 3, 6, 5631,
                                                                       1905, 5649, 163, 169,
                                                                       2121, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6099, 3, 6, 5667,
                                                                       1923, 5703, 181, 191,
                                                                       2139, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6159, 3, 6, 5703,
                                                                       1941, 5739, 191, 201,
                                                                       2169, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6219, 3, 6, 5739,
                                                                       1959, 5775, 201, 211,
                                                                       2199, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6279, 3, 6, 5775,
                                                                       1977, 5811, 211, 221,
                                                                       2229, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6339, 3, 6, 5811,
                                                                       1995, 5847, 221, 231,
                                                                       2259, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6399, 3, 6, 5883,
                                                                       2031, 5919, 251, 261,
                                                                       2289, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6459, 3, 6, 5919,
                                                                       2049, 5955, 261, 271,
                                                                       2319, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6519, 3, 6, 5955,
                                                                       2067, 5991, 271, 281,
                                                                       2349, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6579, 3, 6, 5991,
                                                                       2085, 6027, 281, 291,
                                                                       2379, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 6639, 3, 6, 6027,
                                                                       2103, 6063, 291, 301,
                                                                       2409, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6699, 3, 6, 6099,
                                                                       2139, 6159, 321, 336,
                                                                       2439, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6789, 3, 6, 6159,
                                                                       2169, 6219, 336, 351,
                                                                       2484, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6879, 3, 6, 6219,
                                                                       2199, 6279, 351, 366,
                                                                       2529, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 6969, 3, 6, 6279,
                                                                       2229, 6339, 366, 381,
                                                                       2574, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7059, 3, 6, 6399,
                                                                       2289, 6459, 411, 426,
                                                                       2619, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7149, 3, 6, 6459,
                                                                       2319, 6519, 426, 441,
                                                                       2664, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7239, 3, 6, 6519,
                                                                       2349, 6579, 441, 456,
                                                                       2709, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 7329, 3, 6, 6579,
                                                                       2379, 6639, 456, 471,
                                                                       2754, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7419, 0, 6, 5319,
                                                                       1749, 5325, 2799, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7437, 0, 6, 5325,
                                                                       1752, 5331, 2808, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7455, 0, 6, 5331,
                                                                       1755, 5337, 2817, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7473, 0, 6, 5337,
                                                                       1758, 5343, 2826, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7491, 0, 6, 5343,
                                                                       1761, 5349, 2835, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7509, 0, 6, 5349,
                                                                       1764, 5355, 2844, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7527, 0, 6, 5355,
                                                                       1767, 5361, 2853, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7545, 0, 6, 5367,
                                                                       1773, 5373, 2862, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7563, 0, 6, 5373,
                                                                       1776, 5379, 2871, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7581, 0, 6, 5379,
                                                                       1779, 5385, 2880, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7599, 0, 6, 5385,
                                                                       1782, 5391, 2889, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7617, 0, 6, 5391,
                                                                       1785, 5397, 2898, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7635, 0, 6, 5397,
                                                                       1788, 5403, 2907, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 7653, 0, 6, 5403,
                                                                       1791, 5409, 2916, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7671, 0, 3, 6,
                                                                       5415, 1797, 5433, 7419,
                                                                       2799, 7437, 543, 552,
                                                                       2925, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7725, 0, 3, 6,
                                                                       5433, 1806, 5451, 7437,
                                                                       2808, 7455, 552, 561,
                                                                       2952, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7779, 0, 3, 6,
                                                                       5451, 1815, 5469, 7455,
                                                                       2817, 7473, 561, 570,
                                                                       2979, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7833, 0, 3, 6,
                                                                       5469, 1824, 5487, 7473,
                                                                       2826, 7491, 570, 579,
                                                                       3006, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7887, 0, 3, 6,
                                                                       5487, 1833, 5505, 7491,
                                                                       2835, 7509, 579, 588,
                                                                       3033, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7941, 0, 3, 6,
                                                                       5505, 1842, 5523, 7509,
                                                                       2844, 7527, 588, 597,
                                                                       3060, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 7995, 0, 3, 6,
                                                                       5541, 1860, 5559, 7545,
                                                                       2862, 7563, 615, 624,
                                                                       3087, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8049, 0, 3, 6,
                                                                       5559, 1869, 5577, 7563,
                                                                       2871, 7581, 624, 633,
                                                                       3114, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8103, 0, 3, 6,
                                                                       5577, 1878, 5595, 7581,
                                                                       2880, 7599, 633, 642,
                                                                       3141, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8157, 0, 3, 6,
                                                                       5595, 1887, 5613, 7599,
                                                                       2889, 7617, 642, 651,
                                                                       3168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8211, 0, 3, 6,
                                                                       5613, 1896, 5631, 7617,
                                                                       2898, 7635, 651, 660,
                                                                       3195, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8265, 0, 3, 6,
                                                                       5631, 1905, 5649, 7635,
                                                                       2907, 7653, 660, 669,
                                                                       3222, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 8319, 0, 3, 6,
                                                                       5667, 1923, 5703, 7671,
                                                                       2925, 7725, 687, 705,
                                                                       3249, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 8427, 0, 3, 6,
                                                                       5703, 1941, 5739, 7725,
                                                                       2952, 7779, 705, 723,
                                                                       3303, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 8535, 0, 3, 6,
                                                                       5739, 1959, 5775, 7779,
                                                                       2979, 7833, 723, 741,
                                                                       3357, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 8643, 0, 3, 6,
                                                                       5775, 1977, 5811, 7833,
                                                                       3006, 7887, 741, 759,
                                                                       3411, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 8751, 0, 3, 6,
                                                                       5811, 1995, 5847, 7887,
                                                                       3033, 7941, 759, 777,
                                                                       3465, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 8859, 0, 3, 6,
                                                                       5883, 2031, 5919, 7995,
                                                                       3087, 8049, 813, 831,
                                                                       3519, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 8967, 0, 3, 6,
                                                                       5919, 2049, 5955, 8049,
                                                                       3114, 8103, 831, 849,
                                                                       3573, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 9075, 0, 3, 6,
                                                                       5955, 2067, 5991, 8103,
                                                                       3141, 8157, 849, 867,
                                                                       3627, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 9183, 0, 3, 6,
                                                                       5991, 2085, 6027, 8157,
                                                                       3168, 8211, 867, 885,
                                                                       3681, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 9291, 0, 3, 6,
                                                                       6027, 2103, 6063, 8211,
                                                                       3195, 8265, 885, 903,
                                                                       3735, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 9399, 0, 3, 6,
                                                                       6099, 2139, 6159, 7671,
                                                                       7725, 8319, 3249, 8427,
                                                                       939, 969, 3789, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 9579, 0, 3, 6,
                                                                       6159, 2169, 6219, 7725,
                                                                       7779, 8427, 3303, 8535,
                                                                       969, 999, 3879, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 9759, 0, 3, 6,
                                                                       6219, 2199, 6279, 7779,
                                                                       7833, 8535, 3357, 8643,
                                                                       999, 1029, 3969, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 9939, 0, 3, 6,
                                                                       6279, 2229, 6339, 7833,
                                                                       7887, 8643, 3411, 8751,
                                                                       1029, 1059, 4059, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 10119, 0, 3, 6,
                                                                       6399, 2289, 6459, 7995,
                                                                       8049, 8859, 3519, 8967,
                                                                       1119, 1149, 4149, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 10299, 0, 3, 6,
                                                                       6459, 2319, 6519, 8049,
                                                                       8103, 8967, 3573, 9075,
                                                                       1149, 1179, 4239, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 10479, 0, 3, 6,
                                                                       6519, 2349, 6579, 8103,
                                                                       8157, 9075, 3627, 9183,
                                                                       1179, 1209, 4329, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 10659, 0, 3, 6,
                                                                       6579, 2379, 6639, 8157,
                                                                       8211, 9183, 3681, 9291,
                                                                       1209, 1239, 4419, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 10839, 0, 3, 6,
                                                                       6699, 2439, 6789, 8319,
                                                                       8427, 9399, 3789, 9579,
                                                                       1299, 1344, 4509, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 11109, 0, 3, 6,
                                                                       6789, 2484, 6879, 8427,
                                                                       8535, 9579, 3879, 9759,
                                                                       1344, 1389, 4644, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 11379, 0, 3, 6,
                                                                       6879, 2529, 6969, 8535,
                                                                       8643, 9759, 3969, 9939,
                                                                       1389, 1434, 4779, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 11649, 0, 3, 6,
                                                                       7059, 2619, 7149, 8859,
                                                                       8967, 10119, 4149, 10299,
                                                                       1524, 1569, 4914, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 11919, 0, 3, 6,
                                                                       7149, 2664, 7239, 8967,
                                                                       9075, 10299, 4239, 10479,
                                                                       1569, 1614, 5049, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 12189, 0, 3, 6,
                                                                       7239, 2709, 7329, 9075,
                                                                       9183, 10479, 4329, 10659,
                                                                       1614, 1659, 5184, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12459, 6, 1749,
                                                                       1752, 5331, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12469, 6, 1752,
                                                                       1755, 5337, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12479, 6, 1755,
                                                                       1758, 5343, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12489, 6, 1758,
                                                                       1761, 5349, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12499, 6, 1761,
                                                                       1764, 5355, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12509, 6, 1764,
                                                                       1767, 5361, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12519, 6, 1773,
                                                                       1776, 5379, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12529, 6, 1776,
                                                                       1779, 5385, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12539, 6, 1779,
                                                                       1782, 5391, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12549, 6, 1782,
                                                                       1785, 5397, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12559, 6, 1785,
                                                                       1788, 5403, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 12569, 6, 1788,
                                                                       1791, 5409, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12579, 3, 6,
                                                                       12459, 5331, 12469, 5451,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12609, 3, 6,
                                                                       12469, 5337, 12479, 5469,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12639, 3, 6,
                                                                       12479, 5343, 12489, 5487,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12669, 3, 6,
                                                                       12489, 5349, 12499, 5505,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12699, 3, 6,
                                                                       12499, 5355, 12509, 5523,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12729, 3, 6,
                                                                       12519, 5379, 12529, 5577,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12759, 3, 6,
                                                                       12529, 5385, 12539, 5595,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12789, 3, 6,
                                                                       12539, 5391, 12549, 5613,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12819, 3, 6,
                                                                       12549, 5397, 12559, 5631,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 12849, 3, 6,
                                                                       12559, 5403, 12569, 5649,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12879, 3, 6,
                                                                       12579, 5451, 12609, 1923,
                                                                       1941, 5739, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12939, 3, 6,
                                                                       12609, 5469, 12639, 1941,
                                                                       1959, 5775, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 12999, 3, 6,
                                                                       12639, 5487, 12669, 1959,
                                                                       1977, 5811, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13059, 3, 6,
                                                                       12669, 5505, 12699, 1977,
                                                                       1995, 5847, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13119, 3, 6,
                                                                       12729, 5577, 12759, 2031,
                                                                       2049, 5955, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13179, 3, 6,
                                                                       12759, 5595, 12789, 2049,
                                                                       2067, 5991, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13239, 3, 6,
                                                                       12789, 5613, 12819, 2067,
                                                                       2085, 6027, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 13299, 3, 6,
                                                                       12819, 5631, 12849, 2085,
                                                                       2103, 6063, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13359, 3, 6,
                                                                       12879, 5739, 12939, 2139,
                                                                       2169, 6219, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13459, 3, 6,
                                                                       12939, 5775, 12999, 2169,
                                                                       2199, 6279, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13559, 3, 6,
                                                                       12999, 5811, 13059, 2199,
                                                                       2229, 6339, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13659, 3, 6,
                                                                       13119, 5955, 13179, 2289,
                                                                       2319, 6519, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13759, 3, 6,
                                                                       13179, 5991, 13239, 2319,
                                                                       2349, 6579, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 13859, 3, 6,
                                                                       13239, 6027, 13299, 2349,
                                                                       2379, 6639, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 13959, 3, 6,
                                                                       13359, 6219, 13459, 2439,
                                                                       2484, 6879, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14109, 3, 6,
                                                                       13459, 6279, 13559, 2484,
                                                                       2529, 6969, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14259, 3, 6,
                                                                       13659, 6519, 13759, 2619,
                                                                       2664, 7239, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 14409, 3, 6,
                                                                       13759, 6579, 13859, 2664,
                                                                       2709, 7329, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14559, 0, 6,
                                                                       12459, 5331, 12469, 2799,
                                                                       2808, 7455, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14589, 0, 6,
                                                                       12469, 5337, 12479, 2808,
                                                                       2817, 7473, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14619, 0, 6,
                                                                       12479, 5343, 12489, 2817,
                                                                       2826, 7491, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14649, 0, 6,
                                                                       12489, 5349, 12499, 2826,
                                                                       2835, 7509, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14679, 0, 6,
                                                                       12499, 5355, 12509, 2835,
                                                                       2844, 7527, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14709, 0, 6,
                                                                       12519, 5379, 12529, 2862,
                                                                       2871, 7581, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14739, 0, 6,
                                                                       12529, 5385, 12539, 2871,
                                                                       2880, 7599, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14769, 0, 6,
                                                                       12539, 5391, 12549, 2880,
                                                                       2889, 7617, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14799, 0, 6,
                                                                       12549, 5397, 12559, 2889,
                                                                       2898, 7635, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 14829, 0, 6,
                                                                       12559, 5403, 12569, 2898,
                                                                       2907, 7653, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 14859, 0, 3, 6,
                                                                       12579, 5451, 12609, 14559,
                                                                       7455, 14589, 2925, 2952,
                                                                       7779, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 14949, 0, 3, 6,
                                                                       12609, 5469, 12639, 14589,
                                                                       7473, 14619, 2952, 2979,
                                                                       7833, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15039, 0, 3, 6,
                                                                       12639, 5487, 12669, 14619,
                                                                       7491, 14649, 2979, 3006,
                                                                       7887, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15129, 0, 3, 6,
                                                                       12669, 5505, 12699, 14649,
                                                                       7509, 14679, 3006, 3033,
                                                                       7941, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15219, 0, 3, 6,
                                                                       12729, 5577, 12759, 14709,
                                                                       7581, 14739, 3087, 3114,
                                                                       8103, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15309, 0, 3, 6,
                                                                       12759, 5595, 12789, 14739,
                                                                       7599, 14769, 3114, 3141,
                                                                       8157, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15399, 0, 3, 6,
                                                                       12789, 5613, 12819, 14769,
                                                                       7617, 14799, 3141, 3168,
                                                                       8211, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 15489, 0, 3, 6,
                                                                       12819, 5631, 12849, 14799,
                                                                       7635, 14829, 3168, 3195,
                                                                       8265, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 15579, 0, 3, 6,
                                                                       12879, 5739, 12939, 14859,
                                                                       7779, 14949, 3249, 3303,
                                                                       8535, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 15759, 0, 3, 6,
                                                                       12939, 5775, 12999, 14949,
                                                                       7833, 15039, 3303, 3357,
                                                                       8643, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 15939, 0, 3, 6,
                                                                       12999, 5811, 13059, 15039,
                                                                       7887, 15129, 3357, 3411,
                                                                       8751, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 16119, 0, 3, 6,
                                                                       13119, 5955, 13179, 15219,
                                                                       8103, 15309, 3519, 3573,
                                                                       9075, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 16299, 0, 3, 6,
                                                                       13179, 5991, 13239, 15309,
                                                                       8157, 15399, 3573, 3627,
                                                                       9183, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 16479, 0, 3, 6,
                                                                       13239, 6027, 13299, 15399,
                                                                       8211, 15489, 3627, 3681,
                                                                       9291, ncols, gamma, p,
                                                                       q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 16659, 0, 3, 6,
                                                                       13359, 6219, 13459, 14859,
                                                                       14949, 15579, 8535, 15759,
                                                                       3789, 3879, 9759, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 16959, 0, 3, 6,
                                                                       13459, 6279, 13559, 14949,
                                                                       15039, 15759, 8643, 15939,
                                                                       3879, 3969, 9939, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 17259, 0, 3, 6,
                                                                       13659, 6519, 13759, 15219,
                                                                       15309, 16119, 9075, 16299,
                                                                       4149, 4239, 10479, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 17559, 0, 3, 6,
                                                                       13759, 6579, 13859, 15309,
                                                                       15399, 16299, 9183, 16479,
                                                                       4239, 4329, 10659, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 17859, 0, 3, 6,
                                                                       13959, 6879, 14109, 15579,
                                                                       15759, 16659, 9759, 16959,
                                                                       4509, 4644, 11379, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 18309, 0, 3, 6,
                                                                       14259, 7239, 14409, 16119,
                                                                       16299, 17259, 10479,
                                                                       17559, 4914, 5049, 12189,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18759, 6, 5319,
                                                                       5325, 12459, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18774, 6, 5325,
                                                                       5331, 12469, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18789, 6, 5331,
                                                                       5337, 12479, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18804, 6, 5337,
                                                                       5343, 12489, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18819, 6, 5343,
                                                                       5349, 12499, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18834, 6, 5349,
                                                                       5355, 12509, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18849, 6, 5367,
                                                                       5373, 12519, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18864, 6, 5373,
                                                                       5379, 12529, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18879, 6, 5379,
                                                                       5385, 12539, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18894, 6, 5385,
                                                                       5391, 12549, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18909, 6, 5391,
                                                                       5397, 12559, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 18924, 6, 5397,
                                                                       5403, 12569, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18939, 3, 6,
                                                                       18759, 12459, 18774, 5415,
                                                                       5433, 12579, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 18984, 3, 6,
                                                                       18774, 12469, 18789, 5433,
                                                                       5451, 12609, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19029, 3, 6,
                                                                       18789, 12479, 18804, 5451,
                                                                       5469, 12639, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19074, 3, 6,
                                                                       18804, 12489, 18819, 5469,
                                                                       5487, 12669, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19119, 3, 6,
                                                                       18819, 12499, 18834, 5487,
                                                                       5505, 12699, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19164, 3, 6,
                                                                       18849, 12519, 18864, 5541,
                                                                       5559, 12729, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19209, 3, 6,
                                                                       18864, 12529, 18879, 5559,
                                                                       5577, 12759, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19254, 3, 6,
                                                                       18879, 12539, 18894, 5577,
                                                                       5595, 12789, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19299, 3, 6,
                                                                       18894, 12549, 18909, 5595,
                                                                       5613, 12819, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 19344, 3, 6,
                                                                       18909, 12559, 18924, 5613,
                                                                       5631, 12849, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19389, 3, 6,
                                                                       18939, 12579, 18984, 5667,
                                                                       5703, 12879, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19479, 3, 6,
                                                                       18984, 12609, 19029, 5703,
                                                                       5739, 12939, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19569, 3, 6,
                                                                       19029, 12639, 19074, 5739,
                                                                       5775, 12999, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19659, 3, 6,
                                                                       19074, 12669, 19119, 5775,
                                                                       5811, 13059, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19749, 3, 6,
                                                                       19164, 12729, 19209, 5883,
                                                                       5919, 13119, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19839, 3, 6,
                                                                       19209, 12759, 19254, 5919,
                                                                       5955, 13179, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 19929, 3, 6,
                                                                       19254, 12789, 19299, 5955,
                                                                       5991, 13239, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 20019, 3, 6,
                                                                       19299, 12819, 19344, 5991,
                                                                       6027, 13299, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20109, 3, 6,
                                                                       19389, 12879, 19479, 6099,
                                                                       6159, 13359, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20259, 3, 6,
                                                                       19479, 12939, 19569, 6159,
                                                                       6219, 13459, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20409, 3, 6,
                                                                       19569, 12999, 19659, 6219,
                                                                       6279, 13559, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20559, 3, 6,
                                                                       19749, 13119, 19839, 6399,
                                                                       6459, 13659, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20709, 3, 6,
                                                                       19839, 13179, 19929, 6459,
                                                                       6519, 13759, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 20859, 3, 6,
                                                                       19929, 13239, 20019, 6519,
                                                                       6579, 13859, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 21009, 3, 6,
                                                                       20109, 13359, 20259, 6699,
                                                                       6789, 13959, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 21234, 3, 6,
                                                                       20259, 13459, 20409, 6789,
                                                                       6879, 14109, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 21459, 3, 6,
                                                                       20559, 13659, 20709, 7059,
                                                                       7149, 14259, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 21684, 3, 6,
                                                                       20709, 13759, 20859, 7149,
                                                                       7239, 14409, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 21909, 0, 6,
                                                                       18759, 12459, 18774, 7419,
                                                                       7437, 14559, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 21954, 0, 6,
                                                                       18774, 12469, 18789, 7437,
                                                                       7455, 14589, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 21999, 0, 6,
                                                                       18789, 12479, 18804, 7455,
                                                                       7473, 14619, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 22044, 0, 6,
                                                                       18804, 12489, 18819, 7473,
                                                                       7491, 14649, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 22089, 0, 6,
                                                                       18819, 12499, 18834, 7491,
                                                                       7509, 14679, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 22134, 0, 6,
                                                                       18849, 12519, 18864, 7545,
                                                                       7563, 14709, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 22179, 0, 6,
                                                                       18864, 12529, 18879, 7563,
                                                                       7581, 14739, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 22224, 0, 6,
                                                                       18879, 12539, 18894, 7581,
                                                                       7599, 14769, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 22269, 0, 6,
                                                                       18894, 12549, 18909, 7599,
                                                                       7617, 14799, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 22314, 0, 6,
                                                                       18909, 12559, 18924, 7617,
                                                                       7635, 14829, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 22359, 0, 3, 6,
                                                                       18939, 12579, 18984,
                                                                       21909, 14559, 21954, 7671,
                                                                       7725, 14859, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 22494, 0, 3, 6,
                                                                       18984, 12609, 19029,
                                                                       21954, 14589, 21999, 7725,
                                                                       7779, 14949, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 22629, 0, 3, 6,
                                                                       19029, 12639, 19074,
                                                                       21999, 14619, 22044, 7779,
                                                                       7833, 15039, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 22764, 0, 3, 6,
                                                                       19074, 12669, 19119,
                                                                       22044, 14649, 22089, 7833,
                                                                       7887, 15129, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 22899, 0, 3, 6,
                                                                       19164, 12729, 19209,
                                                                       22134, 14709, 22179, 7995,
                                                                       8049, 15219, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 23034, 0, 3, 6,
                                                                       19209, 12759, 19254,
                                                                       22179, 14739, 22224, 8049,
                                                                       8103, 15309, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 23169, 0, 3, 6,
                                                                       19254, 12789, 19299,
                                                                       22224, 14769, 22269, 8103,
                                                                       8157, 15399, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 23304, 0, 3, 6,
                                                                       19299, 12819, 19344,
                                                                       22269, 14799, 22314, 8157,
                                                                       8211, 15489, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 23439, 0, 3, 6,
                                                                       19389, 12879, 19479,
                                                                       22359, 14859, 22494, 8319,
                                                                       8427, 15579, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 23709, 0, 3, 6,
                                                                       19479, 12939, 19569,
                                                                       22494, 14949, 22629, 8427,
                                                                       8535, 15759, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 23979, 0, 3, 6,
                                                                       19569, 12999, 19659,
                                                                       22629, 15039, 22764, 8535,
                                                                       8643, 15939, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 24249, 0, 3, 6,
                                                                       19749, 13119, 19839,
                                                                       22899, 15219, 23034, 8859,
                                                                       8967, 16119, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 24519, 0, 3, 6,
                                                                       19839, 13179, 19929,
                                                                       23034, 15309, 23169, 8967,
                                                                       9075, 16299, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 24789, 0, 3, 6,
                                                                       19929, 13239, 20019,
                                                                       23169, 15399, 23304, 9075,
                                                                       9183, 16479, ncols, gamma,
                                                                       p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 25059, 0, 3, 6,
                                                                       20109, 13359, 20259,
                                                                       22359, 22494, 23439,
                                                                       15579, 23709, 9399, 9579,
                                                                       16659, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 25509, 0, 3, 6,
                                                                       20259, 13459, 20409,
                                                                       22494, 22629, 23709,
                                                                       15759, 23979, 9579, 9759,
                                                                       16959, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 25959, 0, 3, 6,
                                                                       20559, 13659, 20709,
                                                                       22899, 23034, 24249,
                                                                       16119, 24519, 10119,
                                                                       10299, 17259, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 26409, 0, 3, 6,
                                                                       20709, 13759, 20859,
                                                                       23034, 23169, 24519,
                                                                       16299, 24789, 10299,
                                                                       10479, 17559, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 26859, 0, 3, 6,
                                                                       21009, 13959, 21234,
                                                                       23439, 23709, 25059,
                                                                       16659, 25509, 10839,
                                                                       11109, 17859, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 27534, 0, 3, 6,
                                                                       21459, 14259, 21684,
                                                                       24249, 24519, 25959,
                                                                       17259, 26409, 11649,
                                                                       11919, 18309, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 28209, 27534, 1, 225, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 28434, 27534, 1, 225, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 28659, 27534, 1, 225, ncols, alpha);

                    simdgeo::geom_s_x(buffer, 28884, 26859, 1, 225, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 29109, 26859, 1, 225, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 29334, 26859, 1, 225, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 29559, 28209, 1350, ncols);
                }
            }
        }

        simdtrf::transform_g_inner(buffer, 30909, 29559, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 30909, 9, nmax);

        simdtrf::transform_g_inner(buffer, 30909, 29784, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 81 * nvalues + n * npairs, nvalues, buffer, 30909, 9,
                                   nmax);

        simdtrf::transform_g_inner(buffer, 30909, 30009, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 162 * nvalues + n * npairs, nvalues, buffer, 30909,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 30909, 30234, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 243 * nvalues + n * npairs, nvalues, buffer, 30909,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 30909, 30459, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 324 * nvalues + n * npairs, nvalues, buffer, 30909,
                                   9, nmax);

        simdtrf::transform_g_inner(buffer, 30909, 30684, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 405 * nvalues + n * npairs, nvalues, buffer, 30909,
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
