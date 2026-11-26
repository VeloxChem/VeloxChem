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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecSGH.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGS.hpp"
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
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformH.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_sgh_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_sgh_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 52640, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 594 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 52640, 50585, 1890, dimensions);

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
                                                            4, 5, 6, 7, 8, 9, 10}, ncols, fj,
                                                            i * nprim_b + j, fq, omega);

                    simdfunc::compute_t3c_boys_function(buffer, coordinates, 20, 6, {1, 2, 3, 4,
                                                        5, 6, 7, 8, 9, 10}, ncols, fj,
                                                        i * nprim_b + j, fq);

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

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 501, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 504, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 507, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 510, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 513, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 516, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 519, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 522, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 525, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 528, 0, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 531, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 534, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 537, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 540, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 543, 0, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 546, 0, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 549, 0, 6, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 552, 0, 6, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 555, 0, 6, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 564, 0, 6, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 573, 0, 6, 12, 13,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 582, 0, 6, 13, 14,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 591, 0, 6, 14, 15,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 600, 0, 6, 15, 16,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 609, 0, 6, 16, 17,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 618, 0, 6, 17, 18,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 627, 0, 6, 21, 22,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 636, 0, 6, 22, 23,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 645, 0, 6, 23, 24,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 654, 0, 6, 24, 25,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 663, 0, 6, 25, 26,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 672, 0, 6, 26, 27,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 681, 0, 6, 27, 28,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 690, 0, 6, 28, 29,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 699, 0, 3, 6, 31,
                                                                       34, 85, 91, 555, 564,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 717, 0, 3, 6, 34,
                                                                       37, 91, 97, 564, 573,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 735, 0, 3, 6, 37,
                                                                       40, 97, 103, 573, 582,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 753, 0, 3, 6, 40,
                                                                       43, 103, 109, 582, 591,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 771, 0, 3, 6, 43,
                                                                       46, 109, 115, 591, 600,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 789, 0, 3, 6, 46,
                                                                       49, 115, 121, 600, 609,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 807, 0, 3, 6, 49,
                                                                       52, 121, 127, 609, 618,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 825, 0, 3, 6, 58,
                                                                       61, 133, 139, 627, 636,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 843, 0, 3, 6, 61,
                                                                       64, 139, 145, 636, 645,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 861, 0, 3, 6, 64,
                                                                       67, 145, 151, 645, 654,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 879, 0, 3, 6, 67,
                                                                       70, 151, 157, 654, 663,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 897, 0, 3, 6, 70,
                                                                       73, 157, 163, 663, 672,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 915, 0, 3, 6, 73,
                                                                       76, 163, 169, 672, 681,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 933, 0, 3, 6, 76,
                                                                       79, 169, 175, 681, 690,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 951, 0, 3, 6, 85,
                                                                       91, 181, 191, 699, 717,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 981, 0, 3, 6, 91,
                                                                       97, 191, 201, 717, 735,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1011, 0, 3, 6, 97,
                                                                       103, 201, 211, 735, 753,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1041, 0, 3, 6,
                                                                       103, 109, 211, 221, 753,
                                                                       771, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1071, 0, 3, 6,
                                                                       109, 115, 221, 231, 771,
                                                                       789, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1101, 0, 3, 6,
                                                                       115, 121, 231, 241, 789,
                                                                       807, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1131, 0, 3, 6,
                                                                       133, 139, 251, 261, 825,
                                                                       843, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1161, 0, 3, 6,
                                                                       139, 145, 261, 271, 843,
                                                                       861, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1191, 0, 3, 6,
                                                                       145, 151, 271, 281, 861,
                                                                       879, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1221, 0, 3, 6,
                                                                       151, 157, 281, 291, 879,
                                                                       897, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1251, 0, 3, 6,
                                                                       157, 163, 291, 301, 897,
                                                                       915, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1281, 0, 3, 6,
                                                                       163, 169, 301, 311, 915,
                                                                       933, ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1311, 0, 3, 6,
                                                                       181, 191, 321, 336, 951,
                                                                       981, ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1356, 0, 3, 6,
                                                                       191, 201, 336, 351, 981,
                                                                       1011, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1401, 0, 3, 6,
                                                                       201, 211, 351, 366, 1011,
                                                                       1041, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1446, 0, 3, 6,
                                                                       211, 221, 366, 381, 1041,
                                                                       1071, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1491, 0, 3, 6,
                                                                       221, 231, 381, 396, 1071,
                                                                       1101, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1536, 0, 3, 6,
                                                                       251, 261, 411, 426, 1131,
                                                                       1161, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1581, 0, 3, 6,
                                                                       261, 271, 426, 441, 1161,
                                                                       1191, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1626, 0, 3, 6,
                                                                       271, 281, 441, 456, 1191,
                                                                       1221, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1671, 0, 3, 6,
                                                                       281, 291, 456, 471, 1221,
                                                                       1251, ncols, gamma, p,
                                                                       q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1716, 0, 3, 6,
                                                                       291, 301, 471, 486, 1251,
                                                                       1281, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1761, 6, 10,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1764, 6, 11,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1767, 6, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1770, 6, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1773, 6, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1776, 6, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1779, 6, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1782, 6, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1785, 6, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1788, 6, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1791, 6, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1794, 6, 22,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1797, 6, 23,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1800, 6, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1803, 6, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1806, 6, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1809, 6, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1812, 6, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1815, 6, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1818, 6, 30,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1821, 6, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1830, 6, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1839, 6, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1848, 6, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1857, 6, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1866, 6, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1875, 6, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1884, 6, 23, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1893, 6, 24, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1902, 6, 25, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1911, 6, 26, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1920, 6, 27, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1929, 6, 28, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1938, 6, 29, 82,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1947, 6, 31, 85,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1965, 6, 34, 91,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1983, 6, 37, 97,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2001, 6, 40, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2019, 6, 43, 109,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2037, 6, 46, 115,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2055, 6, 49, 121,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2073, 6, 52, 127,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2091, 6, 58, 133,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2109, 6, 61, 139,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2127, 6, 64, 145,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2145, 6, 67, 151,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2163, 6, 70, 157,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2181, 6, 73, 163,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2199, 6, 76, 169,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 2217, 6, 79, 175,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2235, 6, 85, 181,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2265, 6, 91, 191,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2295, 6, 97, 201,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2325, 6, 103, 211,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2355, 6, 109, 221,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2385, 6, 115, 231,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2415, 6, 121, 241,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2445, 6, 133, 251,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2475, 6, 139, 261,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2505, 6, 145, 271,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2535, 6, 151, 281,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2565, 6, 157, 291,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2595, 6, 163, 301,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2625, 6, 169, 311,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2655, 6, 181, 321,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2700, 6, 191, 336,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2745, 6, 201, 351,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2790, 6, 211, 366,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2835, 6, 221, 381,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2880, 6, 231, 396,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2925, 6, 251, 411,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2970, 6, 261, 426,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3015, 6, 271, 441,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3060, 6, 281, 456,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3105, 6, 291, 471,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 3150, 6, 301, 486,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3195, 6, 10, 501,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3204, 6, 11, 504,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3213, 6, 12, 507,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3222, 6, 13, 510,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3231, 6, 14, 513,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3240, 6, 15, 516,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3249, 6, 16, 519,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3258, 6, 17, 522,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3267, 6, 18, 525,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3276, 6, 21, 528,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3285, 6, 22, 531,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3294, 6, 23, 534,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3303, 6, 24, 537,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3312, 6, 25, 540,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3321, 6, 26, 543,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3330, 6, 27, 546,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3339, 6, 28, 549,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 3348, 6, 29, 552,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3357, 6, 31, 501,
                                                                       555, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3384, 6, 34, 504,
                                                                       564, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3411, 6, 37, 507,
                                                                       573, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3438, 6, 40, 510,
                                                                       582, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3465, 6, 43, 513,
                                                                       591, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3492, 6, 46, 516,
                                                                       600, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3519, 6, 49, 519,
                                                                       609, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3546, 6, 52, 522,
                                                                       618, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3573, 6, 58, 528,
                                                                       627, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3600, 6, 61, 531,
                                                                       636, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3627, 6, 64, 534,
                                                                       645, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3654, 6, 67, 537,
                                                                       654, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3681, 6, 70, 540,
                                                                       663, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3708, 6, 73, 543,
                                                                       672, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3735, 6, 76, 546,
                                                                       681, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3762, 6, 79, 549,
                                                                       690, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3789, 3, 6, 85,
                                                                       3357, 555, 3384, 699,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3843, 3, 6, 91,
                                                                       3384, 564, 3411, 717,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3897, 3, 6, 97,
                                                                       3411, 573, 3438, 735,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3951, 3, 6, 103,
                                                                       3438, 582, 3465, 753,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4005, 3, 6, 109,
                                                                       3465, 591, 3492, 771,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4059, 3, 6, 115,
                                                                       3492, 600, 3519, 789,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4113, 3, 6, 121,
                                                                       3519, 609, 3546, 807,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4167, 3, 6, 133,
                                                                       3573, 627, 3600, 825,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4221, 3, 6, 139,
                                                                       3600, 636, 3627, 843,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4275, 3, 6, 145,
                                                                       3627, 645, 3654, 861,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4329, 3, 6, 151,
                                                                       3654, 654, 3681, 879,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4383, 3, 6, 157,
                                                                       3681, 663, 3708, 897,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4437, 3, 6, 163,
                                                                       3708, 672, 3735, 915,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 4491, 3, 6, 169,
                                                                       3735, 681, 3762, 933,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 4545, 3, 6, 181,
                                                                       3789, 699, 3843, 951,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 4635, 3, 6, 191,
                                                                       3843, 717, 3897, 981,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 4725, 3, 6, 201,
                                                                       3897, 735, 3951, 1011,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 4815, 3, 6, 211,
                                                                       3951, 753, 4005, 1041,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 4905, 3, 6, 221,
                                                                       4005, 771, 4059, 1071,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 4995, 3, 6, 231,
                                                                       4059, 789, 4113, 1101,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5085, 3, 6, 251,
                                                                       4167, 825, 4221, 1131,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5175, 3, 6, 261,
                                                                       4221, 843, 4275, 1161,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5265, 3, 6, 271,
                                                                       4275, 861, 4329, 1191,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5355, 3, 6, 281,
                                                                       4329, 879, 4383, 1221,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5445, 3, 6, 291,
                                                                       4383, 897, 4437, 1251,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 5535, 3, 6, 301,
                                                                       4437, 915, 4491, 1281,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 5625, 3, 6, 321,
                                                                       4545, 951, 4635, 1311,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 5760, 3, 6, 336,
                                                                       4635, 981, 4725, 1356,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 5895, 3, 6, 351,
                                                                       4725, 1011, 4815, 1401,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 6030, 3, 6, 366,
                                                                       4815, 1041, 4905, 1446,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 6165, 3, 6, 381,
                                                                       4905, 1071, 4995, 1491,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 6300, 3, 6, 411,
                                                                       5085, 1131, 5175, 1536,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 6435, 3, 6, 426,
                                                                       5175, 1161, 5265, 1581,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 6570, 3, 6, 441,
                                                                       5265, 1191, 5355, 1626,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 6705, 3, 6, 456,
                                                                       5355, 1221, 5445, 1671,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 6840, 3, 6, 471,
                                                                       5445, 1251, 5535, 1716,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6975, 6, 10, 11,
                                                                       1767, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6981, 6, 11, 12,
                                                                       1770, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6987, 6, 12, 13,
                                                                       1773, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6993, 6, 13, 14,
                                                                       1776, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 6999, 6, 14, 15,
                                                                       1779, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7005, 6, 15, 16,
                                                                       1782, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7011, 6, 16, 17,
                                                                       1785, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7017, 6, 17, 18,
                                                                       1788, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7023, 6, 21, 22,
                                                                       1797, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7029, 6, 22, 23,
                                                                       1800, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7035, 6, 23, 24,
                                                                       1803, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7041, 6, 24, 25,
                                                                       1806, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7047, 6, 25, 26,
                                                                       1809, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7053, 6, 26, 27,
                                                                       1812, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7059, 6, 27, 28,
                                                                       1815, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7065, 6, 28, 29,
                                                                       1818, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7071, 3, 6, 6975,
                                                                       1767, 6981, 1821, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7089, 3, 6, 6981,
                                                                       1770, 6987, 1830, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7107, 3, 6, 6987,
                                                                       1773, 6993, 1839, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7125, 3, 6, 6993,
                                                                       1776, 6999, 1848, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7143, 3, 6, 6999,
                                                                       1779, 7005, 1857, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7161, 3, 6, 7005,
                                                                       1782, 7011, 1866, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7179, 3, 6, 7011,
                                                                       1785, 7017, 1875, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7197, 3, 6, 7023,
                                                                       1797, 7029, 1884, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7215, 3, 6, 7029,
                                                                       1800, 7035, 1893, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7233, 3, 6, 7035,
                                                                       1803, 7041, 1902, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7251, 3, 6, 7041,
                                                                       1806, 7047, 1911, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7269, 3, 6, 7047,
                                                                       1809, 7053, 1920, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7287, 3, 6, 7053,
                                                                       1812, 7059, 1929, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7305, 3, 6, 7059,
                                                                       1815, 7065, 1938, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7323, 3, 6, 7071,
                                                                       1821, 7089, 85, 91, 1983,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7359, 3, 6, 7089,
                                                                       1830, 7107, 91, 97, 2001,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7395, 3, 6, 7107,
                                                                       1839, 7125, 97, 103, 2019,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7431, 3, 6, 7125,
                                                                       1848, 7143, 103, 109,
                                                                       2037, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7467, 3, 6, 7143,
                                                                       1857, 7161, 109, 115,
                                                                       2055, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7503, 3, 6, 7161,
                                                                       1866, 7179, 115, 121,
                                                                       2073, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7539, 3, 6, 7197,
                                                                       1884, 7215, 133, 139,
                                                                       2127, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7575, 3, 6, 7215,
                                                                       1893, 7233, 139, 145,
                                                                       2145, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7611, 3, 6, 7233,
                                                                       1902, 7251, 145, 151,
                                                                       2163, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7647, 3, 6, 7251,
                                                                       1911, 7269, 151, 157,
                                                                       2181, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7683, 3, 6, 7269,
                                                                       1920, 7287, 157, 163,
                                                                       2199, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 7719, 3, 6, 7287,
                                                                       1929, 7305, 163, 169,
                                                                       2217, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7755, 3, 6, 7323,
                                                                       1983, 7359, 181, 191,
                                                                       2295, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7815, 3, 6, 7359,
                                                                       2001, 7395, 191, 201,
                                                                       2325, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7875, 3, 6, 7395,
                                                                       2019, 7431, 201, 211,
                                                                       2355, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7935, 3, 6, 7431,
                                                                       2037, 7467, 211, 221,
                                                                       2385, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 7995, 3, 6, 7467,
                                                                       2055, 7503, 221, 231,
                                                                       2415, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8055, 3, 6, 7539,
                                                                       2127, 7575, 251, 261,
                                                                       2505, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8115, 3, 6, 7575,
                                                                       2145, 7611, 261, 271,
                                                                       2535, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8175, 3, 6, 7611,
                                                                       2163, 7647, 271, 281,
                                                                       2565, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8235, 3, 6, 7647,
                                                                       2181, 7683, 281, 291,
                                                                       2595, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 8295, 3, 6, 7683,
                                                                       2199, 7719, 291, 301,
                                                                       2625, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8355, 3, 6, 7755,
                                                                       2295, 7815, 321, 336,
                                                                       2745, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8445, 3, 6, 7815,
                                                                       2325, 7875, 336, 351,
                                                                       2790, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8535, 3, 6, 7875,
                                                                       2355, 7935, 351, 366,
                                                                       2835, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8625, 3, 6, 7935,
                                                                       2385, 7995, 366, 381,
                                                                       2880, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8715, 3, 6, 8055,
                                                                       2505, 8115, 411, 426,
                                                                       3015, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8805, 3, 6, 8115,
                                                                       2535, 8175, 426, 441,
                                                                       3060, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8895, 3, 6, 8175,
                                                                       2565, 8235, 441, 456,
                                                                       3105, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 8985, 3, 6, 8235,
                                                                       2595, 8295, 456, 471,
                                                                       3150, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9075, 0, 6, 6975,
                                                                       1767, 6981, 3213, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9093, 0, 6, 6981,
                                                                       1770, 6987, 3222, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9111, 0, 6, 6987,
                                                                       1773, 6993, 3231, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9129, 0, 6, 6993,
                                                                       1776, 6999, 3240, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9147, 0, 6, 6999,
                                                                       1779, 7005, 3249, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9165, 0, 6, 7005,
                                                                       1782, 7011, 3258, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9183, 0, 6, 7011,
                                                                       1785, 7017, 3267, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9201, 0, 6, 7023,
                                                                       1797, 7029, 3294, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9219, 0, 6, 7029,
                                                                       1800, 7035, 3303, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9237, 0, 6, 7035,
                                                                       1803, 7041, 3312, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9255, 0, 6, 7041,
                                                                       1806, 7047, 3321, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9273, 0, 6, 7047,
                                                                       1809, 7053, 3330, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9291, 0, 6, 7053,
                                                                       1812, 7059, 3339, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 9309, 0, 6, 7059,
                                                                       1815, 7065, 3348, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9327, 0, 3, 6,
                                                                       7071, 1821, 7089, 9075,
                                                                       3213, 9093, 555, 564,
                                                                       3411, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9381, 0, 3, 6,
                                                                       7089, 1830, 7107, 9093,
                                                                       3222, 9111, 564, 573,
                                                                       3438, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9435, 0, 3, 6,
                                                                       7107, 1839, 7125, 9111,
                                                                       3231, 9129, 573, 582,
                                                                       3465, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9489, 0, 3, 6,
                                                                       7125, 1848, 7143, 9129,
                                                                       3240, 9147, 582, 591,
                                                                       3492, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9543, 0, 3, 6,
                                                                       7143, 1857, 7161, 9147,
                                                                       3249, 9165, 591, 600,
                                                                       3519, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9597, 0, 3, 6,
                                                                       7161, 1866, 7179, 9165,
                                                                       3258, 9183, 600, 609,
                                                                       3546, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9651, 0, 3, 6,
                                                                       7197, 1884, 7215, 9201,
                                                                       3294, 9219, 627, 636,
                                                                       3627, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9705, 0, 3, 6,
                                                                       7215, 1893, 7233, 9219,
                                                                       3303, 9237, 636, 645,
                                                                       3654, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9759, 0, 3, 6,
                                                                       7233, 1902, 7251, 9237,
                                                                       3312, 9255, 645, 654,
                                                                       3681, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9813, 0, 3, 6,
                                                                       7251, 1911, 7269, 9255,
                                                                       3321, 9273, 654, 663,
                                                                       3708, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9867, 0, 3, 6,
                                                                       7269, 1920, 7287, 9273,
                                                                       3330, 9291, 663, 672,
                                                                       3735, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9921, 0, 3, 6,
                                                                       7287, 1929, 7305, 9291,
                                                                       3339, 9309, 672, 681,
                                                                       3762, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 9975, 0, 3, 6,
                                                                       7323, 1983, 7359, 9327,
                                                                       3411, 9381, 699, 717,
                                                                       3897, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 10083, 0, 3, 6,
                                                                       7359, 2001, 7395, 9381,
                                                                       3438, 9435, 717, 735,
                                                                       3951, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 10191, 0, 3, 6,
                                                                       7395, 2019, 7431, 9435,
                                                                       3465, 9489, 735, 753,
                                                                       4005, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 10299, 0, 3, 6,
                                                                       7431, 2037, 7467, 9489,
                                                                       3492, 9543, 753, 771,
                                                                       4059, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 10407, 0, 3, 6,
                                                                       7467, 2055, 7503, 9543,
                                                                       3519, 9597, 771, 789,
                                                                       4113, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 10515, 0, 3, 6,
                                                                       7539, 2127, 7575, 9651,
                                                                       3627, 9705, 825, 843,
                                                                       4275, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 10623, 0, 3, 6,
                                                                       7575, 2145, 7611, 9705,
                                                                       3654, 9759, 843, 861,
                                                                       4329, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 10731, 0, 3, 6,
                                                                       7611, 2163, 7647, 9759,
                                                                       3681, 9813, 861, 879,
                                                                       4383, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 10839, 0, 3, 6,
                                                                       7647, 2181, 7683, 9813,
                                                                       3708, 9867, 879, 897,
                                                                       4437, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 10947, 0, 3, 6,
                                                                       7683, 2199, 7719, 9867,
                                                                       3735, 9921, 897, 915,
                                                                       4491, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 11055, 0, 3, 6,
                                                                       7755, 2295, 7815, 9327,
                                                                       9381, 9975, 3897, 10083,
                                                                       951, 981, 4725, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 11235, 0, 3, 6,
                                                                       7815, 2325, 7875, 9381,
                                                                       9435, 10083, 3951, 10191,
                                                                       981, 1011, 4815, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 11415, 0, 3, 6,
                                                                       7875, 2355, 7935, 9435,
                                                                       9489, 10191, 4005, 10299,
                                                                       1011, 1041, 4905, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 11595, 0, 3, 6,
                                                                       7935, 2385, 7995, 9489,
                                                                       9543, 10299, 4059, 10407,
                                                                       1041, 1071, 4995, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 11775, 0, 3, 6,
                                                                       8055, 2505, 8115, 9651,
                                                                       9705, 10515, 4275, 10623,
                                                                       1131, 1161, 5265, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 11955, 0, 3, 6,
                                                                       8115, 2535, 8175, 9705,
                                                                       9759, 10623, 4329, 10731,
                                                                       1161, 1191, 5355, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 12135, 0, 3, 6,
                                                                       8175, 2565, 8235, 9759,
                                                                       9813, 10731, 4383, 10839,
                                                                       1191, 1221, 5445, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 12315, 0, 3, 6,
                                                                       8235, 2595, 8295, 9813,
                                                                       9867, 10839, 4437, 10947,
                                                                       1221, 1251, 5535, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 12495, 0, 3, 6,
                                                                       8355, 2745, 8445, 9975,
                                                                       10083, 11055, 4725, 11235,
                                                                       1311, 1356, 5895, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 12765, 0, 3, 6,
                                                                       8445, 2790, 8535, 10083,
                                                                       10191, 11235, 4815, 11415,
                                                                       1356, 1401, 6030, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 13035, 0, 3, 6,
                                                                       8535, 2835, 8625, 10191,
                                                                       10299, 11415, 4905, 11595,
                                                                       1401, 1446, 6165, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 13305, 0, 3, 6,
                                                                       8715, 3015, 8805, 10515,
                                                                       10623, 11775, 5265, 11955,
                                                                       1536, 1581, 6570, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 13575, 0, 3, 6,
                                                                       8805, 3060, 8895, 10623,
                                                                       10731, 11955, 5355, 12135,
                                                                       1581, 1626, 6705, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 13845, 0, 3, 6,
                                                                       8895, 3105, 8985, 10731,
                                                                       10839, 12135, 5445, 12315,
                                                                       1626, 1671, 6840, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14115, 6, 1761,
                                                                       1764, 6975, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14125, 6, 1764,
                                                                       1767, 6981, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14135, 6, 1767,
                                                                       1770, 6987, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14145, 6, 1770,
                                                                       1773, 6993, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14155, 6, 1773,
                                                                       1776, 6999, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14165, 6, 1776,
                                                                       1779, 7005, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14175, 6, 1779,
                                                                       1782, 7011, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14185, 6, 1782,
                                                                       1785, 7017, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14195, 6, 1791,
                                                                       1794, 7023, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14205, 6, 1794,
                                                                       1797, 7029, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14215, 6, 1797,
                                                                       1800, 7035, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14225, 6, 1800,
                                                                       1803, 7041, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14235, 6, 1803,
                                                                       1806, 7047, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14245, 6, 1806,
                                                                       1809, 7053, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14255, 6, 1809,
                                                                       1812, 7059, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 14265, 6, 1812,
                                                                       1815, 7065, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14275, 3, 6,
                                                                       14115, 6975, 14125, 7071,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14305, 3, 6,
                                                                       14125, 6981, 14135, 7089,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14335, 3, 6,
                                                                       14135, 6987, 14145, 7107,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14365, 3, 6,
                                                                       14145, 6993, 14155, 7125,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14395, 3, 6,
                                                                       14155, 6999, 14165, 7143,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14425, 3, 6,
                                                                       14165, 7005, 14175, 7161,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14455, 3, 6,
                                                                       14175, 7011, 14185, 7179,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14485, 3, 6,
                                                                       14195, 7023, 14205, 7197,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14515, 3, 6,
                                                                       14205, 7029, 14215, 7215,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14545, 3, 6,
                                                                       14215, 7035, 14225, 7233,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14575, 3, 6,
                                                                       14225, 7041, 14235, 7251,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14605, 3, 6,
                                                                       14235, 7047, 14245, 7269,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14635, 3, 6,
                                                                       14245, 7053, 14255, 7287,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 14665, 3, 6,
                                                                       14255, 7059, 14265, 7305,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14695, 3, 6,
                                                                       14275, 7071, 14305, 1947,
                                                                       1965, 7323, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14755, 3, 6,
                                                                       14305, 7089, 14335, 1965,
                                                                       1983, 7359, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14815, 3, 6,
                                                                       14335, 7107, 14365, 1983,
                                                                       2001, 7395, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14875, 3, 6,
                                                                       14365, 7125, 14395, 2001,
                                                                       2019, 7431, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14935, 3, 6,
                                                                       14395, 7143, 14425, 2019,
                                                                       2037, 7467, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 14995, 3, 6,
                                                                       14425, 7161, 14455, 2037,
                                                                       2055, 7503, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15055, 3, 6,
                                                                       14485, 7197, 14515, 2091,
                                                                       2109, 7539, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15115, 3, 6,
                                                                       14515, 7215, 14545, 2109,
                                                                       2127, 7575, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15175, 3, 6,
                                                                       14545, 7233, 14575, 2127,
                                                                       2145, 7611, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15235, 3, 6,
                                                                       14575, 7251, 14605, 2145,
                                                                       2163, 7647, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15295, 3, 6,
                                                                       14605, 7269, 14635, 2163,
                                                                       2181, 7683, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 15355, 3, 6,
                                                                       14635, 7287, 14665, 2181,
                                                                       2199, 7719, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15415, 3, 6,
                                                                       14695, 7323, 14755, 2235,
                                                                       2265, 7755, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15515, 3, 6,
                                                                       14755, 7359, 14815, 2265,
                                                                       2295, 7815, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15615, 3, 6,
                                                                       14815, 7395, 14875, 2295,
                                                                       2325, 7875, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15715, 3, 6,
                                                                       14875, 7431, 14935, 2325,
                                                                       2355, 7935, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15815, 3, 6,
                                                                       14935, 7467, 14995, 2355,
                                                                       2385, 7995, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 15915, 3, 6,
                                                                       15055, 7539, 15115, 2445,
                                                                       2475, 8055, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 16015, 3, 6,
                                                                       15115, 7575, 15175, 2475,
                                                                       2505, 8115, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 16115, 3, 6,
                                                                       15175, 7611, 15235, 2505,
                                                                       2535, 8175, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 16215, 3, 6,
                                                                       15235, 7647, 15295, 2535,
                                                                       2565, 8235, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 16315, 3, 6,
                                                                       15295, 7683, 15355, 2565,
                                                                       2595, 8295, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16415, 3, 6,
                                                                       15415, 7755, 15515, 2655,
                                                                       2700, 8355, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16565, 3, 6,
                                                                       15515, 7815, 15615, 2700,
                                                                       2745, 8445, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16715, 3, 6,
                                                                       15615, 7875, 15715, 2745,
                                                                       2790, 8535, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 16865, 3, 6,
                                                                       15715, 7935, 15815, 2790,
                                                                       2835, 8625, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17015, 3, 6,
                                                                       15915, 8055, 16015, 2925,
                                                                       2970, 8715, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17165, 3, 6,
                                                                       16015, 8115, 16115, 2970,
                                                                       3015, 8805, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17315, 3, 6,
                                                                       16115, 8175, 16215, 3015,
                                                                       3060, 8895, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 17465, 3, 6,
                                                                       16215, 8235, 16315, 3060,
                                                                       3105, 8985, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17615, 0, 6,
                                                                       14115, 6975, 14125, 3195,
                                                                       3204, 9075, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17645, 0, 6,
                                                                       14125, 6981, 14135, 3204,
                                                                       3213, 9093, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17675, 0, 6,
                                                                       14135, 6987, 14145, 3213,
                                                                       3222, 9111, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17705, 0, 6,
                                                                       14145, 6993, 14155, 3222,
                                                                       3231, 9129, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17735, 0, 6,
                                                                       14155, 6999, 14165, 3231,
                                                                       3240, 9147, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17765, 0, 6,
                                                                       14165, 7005, 14175, 3240,
                                                                       3249, 9165, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17795, 0, 6,
                                                                       14175, 7011, 14185, 3249,
                                                                       3258, 9183, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17825, 0, 6,
                                                                       14195, 7023, 14205, 3276,
                                                                       3285, 9201, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17855, 0, 6,
                                                                       14205, 7029, 14215, 3285,
                                                                       3294, 9219, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17885, 0, 6,
                                                                       14215, 7035, 14225, 3294,
                                                                       3303, 9237, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17915, 0, 6,
                                                                       14225, 7041, 14235, 3303,
                                                                       3312, 9255, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17945, 0, 6,
                                                                       14235, 7047, 14245, 3312,
                                                                       3321, 9273, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 17975, 0, 6,
                                                                       14245, 7053, 14255, 3321,
                                                                       3330, 9291, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18005, 0, 6,
                                                                       14255, 7059, 14265, 3330,
                                                                       3339, 9309, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 18035, 0, 3, 6,
                                                                       14275, 7071, 14305, 17615,
                                                                       9075, 17645, 3357, 3384,
                                                                       9327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 18125, 0, 3, 6,
                                                                       14305, 7089, 14335, 17645,
                                                                       9093, 17675, 3384, 3411,
                                                                       9381, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 18215, 0, 3, 6,
                                                                       14335, 7107, 14365, 17675,
                                                                       9111, 17705, 3411, 3438,
                                                                       9435, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 18305, 0, 3, 6,
                                                                       14365, 7125, 14395, 17705,
                                                                       9129, 17735, 3438, 3465,
                                                                       9489, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 18395, 0, 3, 6,
                                                                       14395, 7143, 14425, 17735,
                                                                       9147, 17765, 3465, 3492,
                                                                       9543, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 18485, 0, 3, 6,
                                                                       14425, 7161, 14455, 17765,
                                                                       9165, 17795, 3492, 3519,
                                                                       9597, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 18575, 0, 3, 6,
                                                                       14485, 7197, 14515, 17825,
                                                                       9201, 17855, 3573, 3600,
                                                                       9651, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 18665, 0, 3, 6,
                                                                       14515, 7215, 14545, 17855,
                                                                       9219, 17885, 3600, 3627,
                                                                       9705, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 18755, 0, 3, 6,
                                                                       14545, 7233, 14575, 17885,
                                                                       9237, 17915, 3627, 3654,
                                                                       9759, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 18845, 0, 3, 6,
                                                                       14575, 7251, 14605, 17915,
                                                                       9255, 17945, 3654, 3681,
                                                                       9813, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 18935, 0, 3, 6,
                                                                       14605, 7269, 14635, 17945,
                                                                       9273, 17975, 3681, 3708,
                                                                       9867, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 19025, 0, 3, 6,
                                                                       14635, 7287, 14665, 17975,
                                                                       9291, 18005, 3708, 3735,
                                                                       9921, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 19115, 0, 3, 6,
                                                                       14695, 7323, 14755, 18035,
                                                                       9327, 18125, 3789, 3843,
                                                                       9975, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 19295, 0, 3, 6,
                                                                       14755, 7359, 14815, 18125,
                                                                       9381, 18215, 3843, 3897,
                                                                       10083, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 19475, 0, 3, 6,
                                                                       14815, 7395, 14875, 18215,
                                                                       9435, 18305, 3897, 3951,
                                                                       10191, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 19655, 0, 3, 6,
                                                                       14875, 7431, 14935, 18305,
                                                                       9489, 18395, 3951, 4005,
                                                                       10299, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 19835, 0, 3, 6,
                                                                       14935, 7467, 14995, 18395,
                                                                       9543, 18485, 4005, 4059,
                                                                       10407, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 20015, 0, 3, 6,
                                                                       15055, 7539, 15115, 18575,
                                                                       9651, 18665, 4167, 4221,
                                                                       10515, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 20195, 0, 3, 6,
                                                                       15115, 7575, 15175, 18665,
                                                                       9705, 18755, 4221, 4275,
                                                                       10623, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 20375, 0, 3, 6,
                                                                       15175, 7611, 15235, 18755,
                                                                       9759, 18845, 4275, 4329,
                                                                       10731, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 20555, 0, 3, 6,
                                                                       15235, 7647, 15295, 18845,
                                                                       9813, 18935, 4329, 4383,
                                                                       10839, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 20735, 0, 3, 6,
                                                                       15295, 7683, 15355, 18935,
                                                                       9867, 19025, 4383, 4437,
                                                                       10947, ncols, gamma, p,
                                                                       q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 20915, 0, 3, 6,
                                                                       15415, 7755, 15515, 18035,
                                                                       18125, 19115, 9975, 19295,
                                                                       4545, 4635, 11055, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 21215, 0, 3, 6,
                                                                       15515, 7815, 15615, 18125,
                                                                       18215, 19295, 10083,
                                                                       19475, 4635, 4725, 11235,
                                                                       ncols, gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 21515, 0, 3, 6,
                                                                       15615, 7875, 15715, 18215,
                                                                       18305, 19475, 10191,
                                                                       19655, 4725, 4815, 11415,
                                                                       ncols, gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 21815, 0, 3, 6,
                                                                       15715, 7935, 15815, 18305,
                                                                       18395, 19655, 10299,
                                                                       19835, 4815, 4905, 11595,
                                                                       ncols, gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 22115, 0, 3, 6,
                                                                       15915, 8055, 16015, 18575,
                                                                       18665, 20015, 10515,
                                                                       20195, 5085, 5175, 11775,
                                                                       ncols, gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 22415, 0, 3, 6,
                                                                       16015, 8115, 16115, 18665,
                                                                       18755, 20195, 10623,
                                                                       20375, 5175, 5265, 11955,
                                                                       ncols, gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 22715, 0, 3, 6,
                                                                       16115, 8175, 16215, 18755,
                                                                       18845, 20375, 10731,
                                                                       20555, 5265, 5355, 12135,
                                                                       ncols, gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 23015, 0, 3, 6,
                                                                       16215, 8235, 16315, 18845,
                                                                       18935, 20555, 10839,
                                                                       20735, 5355, 5445, 12315,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 23315, 0, 3, 6,
                                                                       16415, 8355, 16565, 19115,
                                                                       19295, 20915, 11055,
                                                                       21215, 5625, 5760, 12495,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 23765, 0, 3, 6,
                                                                       16565, 8445, 16715, 19295,
                                                                       19475, 21215, 11235,
                                                                       21515, 5760, 5895, 12765,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 24215, 0, 3, 6,
                                                                       16715, 8535, 16865, 19475,
                                                                       19655, 21515, 11415,
                                                                       21815, 5895, 6030, 13035,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 24665, 0, 3, 6,
                                                                       17015, 8715, 17165, 20015,
                                                                       20195, 22115, 11775,
                                                                       22415, 6300, 6435, 13305,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 25115, 0, 3, 6,
                                                                       17165, 8805, 17315, 20195,
                                                                       20375, 22415, 11955,
                                                                       22715, 6435, 6570, 13575,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 25565, 0, 3, 6,
                                                                       17315, 8895, 17465, 20375,
                                                                       20555, 22715, 12135,
                                                                       23015, 6570, 6705, 13845,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26015, 6, 6975,
                                                                       6981, 14135, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26030, 6, 6981,
                                                                       6987, 14145, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26045, 6, 6987,
                                                                       6993, 14155, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26060, 6, 6993,
                                                                       6999, 14165, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26075, 6, 6999,
                                                                       7005, 14175, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26090, 6, 7005,
                                                                       7011, 14185, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26105, 6, 7023,
                                                                       7029, 14215, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26120, 6, 7029,
                                                                       7035, 14225, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26135, 6, 7035,
                                                                       7041, 14235, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26150, 6, 7041,
                                                                       7047, 14245, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26165, 6, 7047,
                                                                       7053, 14255, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 26180, 6, 7053,
                                                                       7059, 14265, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26195, 3, 6,
                                                                       26015, 14135, 26030, 7071,
                                                                       7089, 14335, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26240, 3, 6,
                                                                       26030, 14145, 26045, 7089,
                                                                       7107, 14365, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26285, 3, 6,
                                                                       26045, 14155, 26060, 7107,
                                                                       7125, 14395, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26330, 3, 6,
                                                                       26060, 14165, 26075, 7125,
                                                                       7143, 14425, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26375, 3, 6,
                                                                       26075, 14175, 26090, 7143,
                                                                       7161, 14455, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26420, 3, 6,
                                                                       26105, 14215, 26120, 7197,
                                                                       7215, 14545, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26465, 3, 6,
                                                                       26120, 14225, 26135, 7215,
                                                                       7233, 14575, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26510, 3, 6,
                                                                       26135, 14235, 26150, 7233,
                                                                       7251, 14605, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26555, 3, 6,
                                                                       26150, 14245, 26165, 7251,
                                                                       7269, 14635, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 26600, 3, 6,
                                                                       26165, 14255, 26180, 7269,
                                                                       7287, 14665, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26645, 3, 6,
                                                                       26195, 14335, 26240, 7323,
                                                                       7359, 14815, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26735, 3, 6,
                                                                       26240, 14365, 26285, 7359,
                                                                       7395, 14875, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26825, 3, 6,
                                                                       26285, 14395, 26330, 7395,
                                                                       7431, 14935, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 26915, 3, 6,
                                                                       26330, 14425, 26375, 7431,
                                                                       7467, 14995, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27005, 3, 6,
                                                                       26420, 14545, 26465, 7539,
                                                                       7575, 15175, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27095, 3, 6,
                                                                       26465, 14575, 26510, 7575,
                                                                       7611, 15235, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27185, 3, 6,
                                                                       26510, 14605, 26555, 7611,
                                                                       7647, 15295, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 27275, 3, 6,
                                                                       26555, 14635, 26600, 7647,
                                                                       7683, 15355, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27365, 3, 6,
                                                                       26645, 14815, 26735, 7755,
                                                                       7815, 15615, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27515, 3, 6,
                                                                       26735, 14875, 26825, 7815,
                                                                       7875, 15715, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27665, 3, 6,
                                                                       26825, 14935, 26915, 7875,
                                                                       7935, 15815, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27815, 3, 6,
                                                                       27005, 15175, 27095, 8055,
                                                                       8115, 16115, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 27965, 3, 6,
                                                                       27095, 15235, 27185, 8115,
                                                                       8175, 16215, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 28115, 3, 6,
                                                                       27185, 15295, 27275, 8175,
                                                                       8235, 16315, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28265, 3, 6,
                                                                       27365, 15615, 27515, 8355,
                                                                       8445, 16715, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28490, 3, 6,
                                                                       27515, 15715, 27665, 8445,
                                                                       8535, 16865, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28715, 3, 6,
                                                                       27815, 16115, 27965, 8715,
                                                                       8805, 17315, ncols, gamma,
                                                                       p, q);

                    compute_prim_sgg_three_center_electron_repulsion_0(buffer, 28940, 3, 6,
                                                                       27965, 16215, 28115, 8805,
                                                                       8895, 17465, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29165, 0, 6,
                                                                       26015, 14135, 26030, 9075,
                                                                       9093, 17675, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29210, 0, 6,
                                                                       26030, 14145, 26045, 9093,
                                                                       9111, 17705, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29255, 0, 6,
                                                                       26045, 14155, 26060, 9111,
                                                                       9129, 17735, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29300, 0, 6,
                                                                       26060, 14165, 26075, 9129,
                                                                       9147, 17765, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29345, 0, 6,
                                                                       26075, 14175, 26090, 9147,
                                                                       9165, 17795, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29390, 0, 6,
                                                                       26105, 14215, 26120, 9201,
                                                                       9219, 17885, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29435, 0, 6,
                                                                       26120, 14225, 26135, 9219,
                                                                       9237, 17915, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29480, 0, 6,
                                                                       26135, 14235, 26150, 9237,
                                                                       9255, 17945, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29525, 0, 6,
                                                                       26150, 14245, 26165, 9255,
                                                                       9273, 17975, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 29570, 0, 6,
                                                                       26165, 14255, 26180, 9273,
                                                                       9291, 18005, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 29615, 0, 3, 6,
                                                                       26195, 14335, 26240,
                                                                       29165, 17675, 29210, 9327,
                                                                       9381, 18215, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 29750, 0, 3, 6,
                                                                       26240, 14365, 26285,
                                                                       29210, 17705, 29255, 9381,
                                                                       9435, 18305, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 29885, 0, 3, 6,
                                                                       26285, 14395, 26330,
                                                                       29255, 17735, 29300, 9435,
                                                                       9489, 18395, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 30020, 0, 3, 6,
                                                                       26330, 14425, 26375,
                                                                       29300, 17765, 29345, 9489,
                                                                       9543, 18485, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 30155, 0, 3, 6,
                                                                       26420, 14545, 26465,
                                                                       29390, 17885, 29435, 9651,
                                                                       9705, 18755, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 30290, 0, 3, 6,
                                                                       26465, 14575, 26510,
                                                                       29435, 17915, 29480, 9705,
                                                                       9759, 18845, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 30425, 0, 3, 6,
                                                                       26510, 14605, 26555,
                                                                       29480, 17945, 29525, 9759,
                                                                       9813, 18935, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 30560, 0, 3, 6,
                                                                       26555, 14635, 26600,
                                                                       29525, 17975, 29570, 9813,
                                                                       9867, 19025, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 30695, 0, 3, 6,
                                                                       26645, 14815, 26735,
                                                                       29615, 18215, 29750, 9975,
                                                                       10083, 19475, ncols,
                                                                       gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 30965, 0, 3, 6,
                                                                       26735, 14875, 26825,
                                                                       29750, 18305, 29885,
                                                                       10083, 10191, 19655,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 31235, 0, 3, 6,
                                                                       26825, 14935, 26915,
                                                                       29885, 18395, 30020,
                                                                       10191, 10299, 19835,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 31505, 0, 3, 6,
                                                                       27005, 15175, 27095,
                                                                       30155, 18755, 30290,
                                                                       10515, 10623, 20375,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 31775, 0, 3, 6,
                                                                       27095, 15235, 27185,
                                                                       30290, 18845, 30425,
                                                                       10623, 10731, 20555,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 32045, 0, 3, 6,
                                                                       27185, 15295, 27275,
                                                                       30425, 18935, 30560,
                                                                       10731, 10839, 20735,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 32315, 0, 3, 6,
                                                                       27365, 15615, 27515,
                                                                       29615, 29750, 30695,
                                                                       19475, 30965, 11055,
                                                                       11235, 21515, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 32765, 0, 3, 6,
                                                                       27515, 15715, 27665,
                                                                       29750, 29885, 30965,
                                                                       19655, 31235, 11235,
                                                                       11415, 21815, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 33215, 0, 3, 6,
                                                                       27815, 16115, 27965,
                                                                       30155, 30290, 31505,
                                                                       20375, 31775, 11775,
                                                                       11955, 22715, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 33665, 0, 3, 6,
                                                                       27965, 16215, 28115,
                                                                       30290, 30425, 31775,
                                                                       20555, 32045, 11955,
                                                                       12135, 23015, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 34115, 0, 3, 6,
                                                                       28265, 16715, 28490,
                                                                       30695, 30965, 32315,
                                                                       21515, 32765, 12495,
                                                                       12765, 24215, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgg_three_center_electron_repulsion_0(buffer, 34790, 0, 3, 6,
                                                                       28715, 17315, 28940,
                                                                       31505, 31775, 33215,
                                                                       22715, 33665, 13305,
                                                                       13575, 25565, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35465, 6, 14115,
                                                                       14125, 26015, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35486, 6, 14125,
                                                                       14135, 26030, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35507, 6, 14135,
                                                                       14145, 26045, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35528, 6, 14145,
                                                                       14155, 26060, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35549, 6, 14155,
                                                                       14165, 26075, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35570, 6, 14165,
                                                                       14175, 26090, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35591, 6, 14195,
                                                                       14205, 26105, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35612, 6, 14205,
                                                                       14215, 26120, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35633, 6, 14215,
                                                                       14225, 26135, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35654, 6, 14225,
                                                                       14235, 26150, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35675, 6, 14235,
                                                                       14245, 26165, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 35696, 6, 14245,
                                                                       14255, 26180, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35717, 3, 6,
                                                                       35465, 26015, 35486,
                                                                       14275, 14305, 26195,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35780, 3, 6,
                                                                       35486, 26030, 35507,
                                                                       14305, 14335, 26240,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35843, 3, 6,
                                                                       35507, 26045, 35528,
                                                                       14335, 14365, 26285,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35906, 3, 6,
                                                                       35528, 26060, 35549,
                                                                       14365, 14395, 26330,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 35969, 3, 6,
                                                                       35549, 26075, 35570,
                                                                       14395, 14425, 26375,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36032, 3, 6,
                                                                       35591, 26105, 35612,
                                                                       14485, 14515, 26420,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36095, 3, 6,
                                                                       35612, 26120, 35633,
                                                                       14515, 14545, 26465,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36158, 3, 6,
                                                                       35633, 26135, 35654,
                                                                       14545, 14575, 26510,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36221, 3, 6,
                                                                       35654, 26150, 35675,
                                                                       14575, 14605, 26555,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 36284, 3, 6,
                                                                       35675, 26165, 35696,
                                                                       14605, 14635, 26600,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36347, 3, 6,
                                                                       35717, 26195, 35780,
                                                                       14695, 14755, 26645,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36473, 3, 6,
                                                                       35780, 26240, 35843,
                                                                       14755, 14815, 26735,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36599, 3, 6,
                                                                       35843, 26285, 35906,
                                                                       14815, 14875, 26825,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36725, 3, 6,
                                                                       35906, 26330, 35969,
                                                                       14875, 14935, 26915,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36851, 3, 6,
                                                                       36032, 26420, 36095,
                                                                       15055, 15115, 27005,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 36977, 3, 6,
                                                                       36095, 26465, 36158,
                                                                       15115, 15175, 27095,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37103, 3, 6,
                                                                       36158, 26510, 36221,
                                                                       15175, 15235, 27185,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 37229, 3, 6,
                                                                       36221, 26555, 36284,
                                                                       15235, 15295, 27275,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37355, 3, 6,
                                                                       36347, 26645, 36473,
                                                                       15415, 15515, 27365,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37565, 3, 6,
                                                                       36473, 26735, 36599,
                                                                       15515, 15615, 27515,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37775, 3, 6,
                                                                       36599, 26825, 36725,
                                                                       15615, 15715, 27665,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 37985, 3, 6,
                                                                       36851, 27005, 36977,
                                                                       15915, 16015, 27815,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38195, 3, 6,
                                                                       36977, 27095, 37103,
                                                                       16015, 16115, 27965,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 38405, 3, 6,
                                                                       37103, 27185, 37229,
                                                                       16115, 16215, 28115,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 38615, 3, 6,
                                                                       37355, 27365, 37565,
                                                                       16415, 16565, 28265,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 38930, 3, 6,
                                                                       37565, 27515, 37775,
                                                                       16565, 16715, 28490,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 39245, 3, 6,
                                                                       37985, 27815, 38195,
                                                                       17015, 17165, 28715,
                                                                       ncols, gamma, p, q);

                    compute_prim_sgh_three_center_electron_repulsion_0(buffer, 39560, 3, 6,
                                                                       38195, 27965, 38405,
                                                                       17165, 17315, 28940,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 39875, 0, 6,
                                                                       35465, 26015, 35486,
                                                                       17615, 17645, 29165,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 39938, 0, 6,
                                                                       35486, 26030, 35507,
                                                                       17645, 17675, 29210,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 40001, 0, 6,
                                                                       35507, 26045, 35528,
                                                                       17675, 17705, 29255,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 40064, 0, 6,
                                                                       35528, 26060, 35549,
                                                                       17705, 17735, 29300,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 40127, 0, 6,
                                                                       35549, 26075, 35570,
                                                                       17735, 17765, 29345,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 40190, 0, 6,
                                                                       35591, 26105, 35612,
                                                                       17825, 17855, 29390,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 40253, 0, 6,
                                                                       35612, 26120, 35633,
                                                                       17855, 17885, 29435,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 40316, 0, 6,
                                                                       35633, 26135, 35654,
                                                                       17885, 17915, 29480,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 40379, 0, 6,
                                                                       35654, 26150, 35675,
                                                                       17915, 17945, 29525,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 40442, 0, 6,
                                                                       35675, 26165, 35696,
                                                                       17945, 17975, 29570,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 40505, 0, 3, 6,
                                                                       35717, 26195, 35780,
                                                                       39875, 29165, 39938,
                                                                       18035, 18125, 29615,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 40694, 0, 3, 6,
                                                                       35780, 26240, 35843,
                                                                       39938, 29210, 40001,
                                                                       18125, 18215, 29750,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 40883, 0, 3, 6,
                                                                       35843, 26285, 35906,
                                                                       40001, 29255, 40064,
                                                                       18215, 18305, 29885,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 41072, 0, 3, 6,
                                                                       35906, 26330, 35969,
                                                                       40064, 29300, 40127,
                                                                       18305, 18395, 30020,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 41261, 0, 3, 6,
                                                                       36032, 26420, 36095,
                                                                       40190, 29390, 40253,
                                                                       18575, 18665, 30155,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 41450, 0, 3, 6,
                                                                       36095, 26465, 36158,
                                                                       40253, 29435, 40316,
                                                                       18665, 18755, 30290,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 41639, 0, 3, 6,
                                                                       36158, 26510, 36221,
                                                                       40316, 29480, 40379,
                                                                       18755, 18845, 30425,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 41828, 0, 3, 6,
                                                                       36221, 26555, 36284,
                                                                       40379, 29525, 40442,
                                                                       18845, 18935, 30560,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 42017, 0, 3, 6,
                                                                       36347, 26645, 36473,
                                                                       40505, 29615, 40694,
                                                                       19115, 19295, 30695,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 42395, 0, 3, 6,
                                                                       36473, 26735, 36599,
                                                                       40694, 29750, 40883,
                                                                       19295, 19475, 30965,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 42773, 0, 3, 6,
                                                                       36599, 26825, 36725,
                                                                       40883, 29885, 41072,
                                                                       19475, 19655, 31235,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 43151, 0, 3, 6,
                                                                       36851, 27005, 36977,
                                                                       41261, 30155, 41450,
                                                                       20015, 20195, 31505,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 43529, 0, 3, 6,
                                                                       36977, 27095, 37103,
                                                                       41450, 30290, 41639,
                                                                       20195, 20375, 31775,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 43907, 0, 3, 6,
                                                                       37103, 27185, 37229,
                                                                       41639, 30425, 41828,
                                                                       20375, 20555, 32045,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 44285, 0, 3, 6,
                                                                       37355, 27365, 37565,
                                                                       40505, 40694, 42017,
                                                                       30695, 42395, 20915,
                                                                       21215, 32315, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 44915, 0, 3, 6,
                                                                       37565, 27515, 37775,
                                                                       40694, 40883, 42395,
                                                                       30965, 42773, 21215,
                                                                       21515, 32765, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 45545, 0, 3, 6,
                                                                       37985, 27815, 38195,
                                                                       41261, 41450, 43151,
                                                                       31505, 43529, 22115,
                                                                       22415, 33215, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 46175, 0, 3, 6,
                                                                       38195, 27965, 38405,
                                                                       41450, 41639, 43529,
                                                                       31775, 43907, 22415,
                                                                       22715, 33665, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgh_three_center_electron_repulsion_0(buffer, 46805, 0, 3, 6,
                                                                       38615, 28265, 38930,
                                                                       42017, 42395, 44285,
                                                                       32315, 44915, 23315,
                                                                       23765, 34115, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgh_three_center_electron_repulsion_0(buffer, 47750, 0, 3, 6,
                                                                       39245, 28715, 39560,
                                                                       43151, 43529, 45545,
                                                                       33215, 46175, 24665,
                                                                       25115, 34790, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 48695, 47750, 1, 315, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 49010, 47750, 1, 315, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 49325, 47750, 1, 315, ncols, alpha);

                    simdgeo::geom_s_x(buffer, 49640, 46805, 1, 315, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 49955, 46805, 1, 315, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 50270, 46805, 1, 315, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 50585, 48695, 1890, ncols);
                }
            }
        }

        simdtrf::transform_h_inner(buffer, 52475, 50585, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 52475, 11, nmax);

        simdtrf::transform_h_inner(buffer, 52475, 50900, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 99 * nvalues + n * npairs, nvalues, buffer, 52475,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 52475, 51215, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 198 * nvalues + n * npairs, nvalues, buffer, 52475,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 52475, 51530, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 297 * nvalues + n * npairs, nvalues, buffer, 52475,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 52475, 51845, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 396 * nvalues + n * npairs, nvalues, buffer, 52475,
                                   11, nmax);

        simdtrf::transform_h_inner(buffer, 52475, 52160, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 495 * nvalues + n * npairs, nvalues, buffer, 52475,
                                   11, nmax);
    }

    for (size_t m = 0; m < 594; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
