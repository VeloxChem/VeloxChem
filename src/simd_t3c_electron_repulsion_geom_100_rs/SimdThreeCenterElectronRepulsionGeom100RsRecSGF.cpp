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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecSGF.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPGS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSGS.hpp"
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
compute_rs_geom_100_sgf_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_sgf_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 69, 3, 6, 10, 11,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 75, 3, 6, 11, 12,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 81, 3, 6, 12, 13,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 87, 3, 6, 13, 14,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 93, 3, 6, 14, 15,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 99, 3, 6, 15, 16,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 105, 3, 6, 19, 20,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 111, 3, 6, 20, 21,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 117, 3, 6, 21, 22,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 123, 3, 6, 22, 23,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 129, 3, 6, 23, 24,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 135, 3, 6, 24, 25,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 141, 3, 6, 27, 30,
                                                                       69, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 151, 3, 6, 30, 33,
                                                                       75, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 161, 3, 6, 33, 36,
                                                                       81, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 171, 3, 6, 36, 39,
                                                                       87, 93, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 181, 3, 6, 39, 42,
                                                                       93, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 191, 3, 6, 48, 51,
                                                                       105, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 201, 3, 6, 51, 54,
                                                                       111, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 211, 3, 6, 54, 57,
                                                                       117, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 221, 3, 6, 57, 60,
                                                                       123, 129, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 231, 3, 6, 60, 63,
                                                                       129, 135, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 241, 3, 6, 69, 75,
                                                                       141, 151, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 256, 3, 6, 75, 81,
                                                                       151, 161, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 271, 3, 6, 81, 87,
                                                                       161, 171, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 286, 3, 6, 87, 93,
                                                                       171, 181, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 301, 3, 6, 105,
                                                                       111, 191, 201, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 316, 3, 6, 111,
                                                                       117, 201, 211, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 331, 3, 6, 117,
                                                                       123, 211, 221, ncols,
                                                                       gamma, p, q);

                    compute_prim_sgs_three_center_electron_repulsion_0(buffer, 346, 3, 6, 123,
                                                                       129, 221, 231, ncols,
                                                                       gamma, p, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 361, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 364, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 367, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 370, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 373, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 376, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 379, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 382, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 385, 0, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 388, 0, 6, 21, 22,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 391, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 394, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 397, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 400, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 403, 0, 6, 10, 11,
                                                                       27, 30, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 412, 0, 6, 11, 12,
                                                                       30, 33, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 421, 0, 6, 12, 13,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 430, 0, 6, 13, 14,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 439, 0, 6, 14, 15,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 448, 0, 6, 15, 16,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 457, 0, 6, 19, 20,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 466, 0, 6, 20, 21,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 475, 0, 6, 21, 22,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 484, 0, 6, 22, 23,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 493, 0, 6, 23, 24,
                                                                       60, 63, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 502, 0, 6, 24, 25,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 511, 0, 3, 6, 27,
                                                                       30, 69, 75, 403, 412,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 529, 0, 3, 6, 30,
                                                                       33, 75, 81, 412, 421,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 547, 0, 3, 6, 33,
                                                                       36, 81, 87, 421, 430,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 565, 0, 3, 6, 36,
                                                                       39, 87, 93, 430, 439,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 583, 0, 3, 6, 39,
                                                                       42, 93, 99, 439, 448,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 601, 0, 3, 6, 48,
                                                                       51, 105, 111, 457, 466,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 619, 0, 3, 6, 51,
                                                                       54, 111, 117, 466, 475,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 637, 0, 3, 6, 54,
                                                                       57, 117, 123, 475, 484,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 655, 0, 3, 6, 57,
                                                                       60, 123, 129, 484, 493,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 673, 0, 3, 6, 60,
                                                                       63, 129, 135, 493, 502,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 691, 0, 3, 6, 69,
                                                                       75, 141, 151, 511, 529,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 721, 0, 3, 6, 75,
                                                                       81, 151, 161, 529, 547,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 751, 0, 3, 6, 81,
                                                                       87, 161, 171, 547, 565,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 781, 0, 3, 6, 87,
                                                                       93, 171, 181, 565, 583,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 811, 0, 3, 6, 105,
                                                                       111, 191, 201, 601, 619,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 841, 0, 3, 6, 111,
                                                                       117, 201, 211, 619, 637,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 871, 0, 3, 6, 117,
                                                                       123, 211, 221, 637, 655,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 901, 0, 3, 6, 123,
                                                                       129, 221, 231, 655, 673,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 931, 0, 3, 6, 141,
                                                                       151, 241, 256, 691, 721,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 976, 0, 3, 6, 151,
                                                                       161, 256, 271, 721, 751,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1021, 0, 3, 6,
                                                                       161, 171, 271, 286, 751,
                                                                       781, ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1066, 0, 3, 6,
                                                                       191, 201, 301, 316, 811,
                                                                       841, ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1111, 0, 3, 6,
                                                                       201, 211, 316, 331, 841,
                                                                       871, ncols, gamma, p, q);

                    compute_prim_pgs_three_center_electron_repulsion_0(buffer, 1156, 0, 3, 6,
                                                                       211, 221, 331, 346, 871,
                                                                       901, ncols, gamma, p, q);

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

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1339, 6, 27, 69,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1357, 6, 30, 75,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1375, 6, 33, 81,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1393, 6, 36, 87,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1411, 6, 39, 93,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1429, 6, 42, 99,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1447, 6, 48, 105,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1465, 6, 51, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1483, 6, 54, 117,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1501, 6, 57, 123,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1519, 6, 60, 129,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1537, 6, 63, 135,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1555, 6, 69, 141,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1585, 6, 75, 151,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1615, 6, 81, 161,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1645, 6, 87, 171,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1675, 6, 93, 181,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1705, 6, 105, 191,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1735, 6, 111, 201,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1765, 6, 117, 211,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1795, 6, 123, 221,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1825, 6, 129, 231,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1855, 6, 141, 241,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1900, 6, 151, 256,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1945, 6, 161, 271,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 1990, 6, 171, 286,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2035, 6, 191, 301,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2080, 6, 201, 316,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2125, 6, 211, 331,
                                                                       ncols, p, q);

                    compute_prim_sgp_three_center_electron_repulsion_0(buffer, 2170, 6, 221, 346,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2215, 6, 10, 361,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2224, 6, 11, 364,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2233, 6, 12, 367,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2242, 6, 13, 370,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2251, 6, 14, 373,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2260, 6, 15, 376,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2269, 6, 16, 379,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2278, 6, 19, 382,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2287, 6, 20, 385,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2296, 6, 21, 388,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2305, 6, 22, 391,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2314, 6, 23, 394,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2323, 6, 24, 397,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2332, 6, 25, 400,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2341, 6, 27, 361,
                                                                       403, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2368, 6, 30, 364,
                                                                       412, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2395, 6, 33, 367,
                                                                       421, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2422, 6, 36, 370,
                                                                       430, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2449, 6, 39, 373,
                                                                       439, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2476, 6, 42, 376,
                                                                       448, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2503, 6, 48, 382,
                                                                       457, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2530, 6, 51, 385,
                                                                       466, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2557, 6, 54, 388,
                                                                       475, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2584, 6, 57, 391,
                                                                       484, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2611, 6, 60, 394,
                                                                       493, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2638, 6, 63, 397,
                                                                       502, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2665, 3, 6, 69,
                                                                       2341, 403, 2368, 511,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2719, 3, 6, 75,
                                                                       2368, 412, 2395, 529,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2773, 3, 6, 81,
                                                                       2395, 421, 2422, 547,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2827, 3, 6, 87,
                                                                       2422, 430, 2449, 565,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2881, 3, 6, 93,
                                                                       2449, 439, 2476, 583,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2935, 3, 6, 105,
                                                                       2503, 457, 2530, 601,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2989, 3, 6, 111,
                                                                       2530, 466, 2557, 619,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3043, 3, 6, 117,
                                                                       2557, 475, 2584, 637,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3097, 3, 6, 123,
                                                                       2584, 484, 2611, 655,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3151, 3, 6, 129,
                                                                       2611, 493, 2638, 673,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3205, 3, 6, 141,
                                                                       2665, 511, 2719, 691,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3295, 3, 6, 151,
                                                                       2719, 529, 2773, 721,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3385, 3, 6, 161,
                                                                       2773, 547, 2827, 751,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3475, 3, 6, 171,
                                                                       2827, 565, 2881, 781,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3565, 3, 6, 191,
                                                                       2935, 601, 2989, 811,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3655, 3, 6, 201,
                                                                       2989, 619, 3043, 841,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3745, 3, 6, 211,
                                                                       3043, 637, 3097, 871,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3835, 3, 6, 221,
                                                                       3097, 655, 3151, 901,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 3925, 3, 6, 241,
                                                                       3205, 691, 3295, 931,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 4060, 3, 6, 256,
                                                                       3295, 721, 3385, 976,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 4195, 3, 6, 271,
                                                                       3385, 751, 3475, 1021,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 4330, 3, 6, 301,
                                                                       3565, 811, 3655, 1066,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 4465, 3, 6, 316,
                                                                       3655, 841, 3745, 1111,
                                                                       ncols, gamma, p, q);

                    compute_prim_pgp_three_center_electron_repulsion_0(buffer, 4600, 3, 6, 331,
                                                                       3745, 871, 3835, 1156,
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

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4987, 3, 6, 4807,
                                                                       1249, 4825, 69, 75, 1375,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5023, 3, 6, 4825,
                                                                       1258, 4843, 75, 81, 1393,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5059, 3, 6, 4843,
                                                                       1267, 4861, 81, 87, 1411,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5095, 3, 6, 4861,
                                                                       1276, 4879, 87, 93, 1429,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5131, 3, 6, 4897,
                                                                       1294, 4915, 105, 111,
                                                                       1483, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5167, 3, 6, 4915,
                                                                       1303, 4933, 111, 117,
                                                                       1501, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5203, 3, 6, 4933,
                                                                       1312, 4951, 117, 123,
                                                                       1519, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5239, 3, 6, 4951,
                                                                       1321, 4969, 123, 129,
                                                                       1537, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5275, 3, 6, 4987,
                                                                       1375, 5023, 141, 151,
                                                                       1615, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5335, 3, 6, 5023,
                                                                       1393, 5059, 151, 161,
                                                                       1645, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5395, 3, 6, 5059,
                                                                       1411, 5095, 161, 171,
                                                                       1675, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5455, 3, 6, 5131,
                                                                       1483, 5167, 191, 201,
                                                                       1765, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5515, 3, 6, 5167,
                                                                       1501, 5203, 201, 211,
                                                                       1795, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5575, 3, 6, 5203,
                                                                       1519, 5239, 211, 221,
                                                                       1825, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5635, 3, 6, 5275,
                                                                       1615, 5335, 241, 256,
                                                                       1945, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5725, 3, 6, 5335,
                                                                       1645, 5395, 256, 271,
                                                                       1990, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5815, 3, 6, 5455,
                                                                       1765, 5515, 301, 316,
                                                                       2125, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgd_three_center_electron_repulsion_0(buffer, 5905, 3, 6, 5515,
                                                                       1795, 5575, 316, 331,
                                                                       2170, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5995, 0, 6, 4735,
                                                                       1207, 4741, 2233, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6013, 0, 6, 4741,
                                                                       1210, 4747, 2242, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6031, 0, 6, 4747,
                                                                       1213, 4753, 2251, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6049, 0, 6, 4753,
                                                                       1216, 4759, 2260, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6067, 0, 6, 4759,
                                                                       1219, 4765, 2269, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6085, 0, 6, 4771,
                                                                       1231, 4777, 2296, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6103, 0, 6, 4777,
                                                                       1234, 4783, 2305, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6121, 0, 6, 4783,
                                                                       1237, 4789, 2314, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6139, 0, 6, 4789,
                                                                       1240, 4795, 2323, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6157, 0, 6, 4795,
                                                                       1243, 4801, 2332, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6175, 0, 3, 6,
                                                                       4807, 1249, 4825, 5995,
                                                                       2233, 6013, 403, 412,
                                                                       2395, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6229, 0, 3, 6,
                                                                       4825, 1258, 4843, 6013,
                                                                       2242, 6031, 412, 421,
                                                                       2422, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6283, 0, 3, 6,
                                                                       4843, 1267, 4861, 6031,
                                                                       2251, 6049, 421, 430,
                                                                       2449, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6337, 0, 3, 6,
                                                                       4861, 1276, 4879, 6049,
                                                                       2260, 6067, 430, 439,
                                                                       2476, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6391, 0, 3, 6,
                                                                       4897, 1294, 4915, 6085,
                                                                       2296, 6103, 457, 466,
                                                                       2557, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6445, 0, 3, 6,
                                                                       4915, 1303, 4933, 6103,
                                                                       2305, 6121, 466, 475,
                                                                       2584, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6499, 0, 3, 6,
                                                                       4933, 1312, 4951, 6121,
                                                                       2314, 6139, 475, 484,
                                                                       2611, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6553, 0, 3, 6,
                                                                       4951, 1321, 4969, 6139,
                                                                       2323, 6157, 484, 493,
                                                                       2638, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6607, 0, 3, 6,
                                                                       4987, 1375, 5023, 6175,
                                                                       2395, 6229, 511, 529,
                                                                       2773, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6715, 0, 3, 6,
                                                                       5023, 1393, 5059, 6229,
                                                                       2422, 6283, 529, 547,
                                                                       2827, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6823, 0, 3, 6,
                                                                       5059, 1411, 5095, 6283,
                                                                       2449, 6337, 547, 565,
                                                                       2881, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6931, 0, 3, 6,
                                                                       5131, 1483, 5167, 6391,
                                                                       2557, 6445, 601, 619,
                                                                       3043, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 7039, 0, 3, 6,
                                                                       5167, 1501, 5203, 6445,
                                                                       2584, 6499, 619, 637,
                                                                       3097, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 7147, 0, 3, 6,
                                                                       5203, 1519, 5239, 6499,
                                                                       2611, 6553, 637, 655,
                                                                       3151, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 7255, 0, 3, 6,
                                                                       5275, 1615, 5335, 6175,
                                                                       6229, 6607, 2773, 6715,
                                                                       691, 721, 3385, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 7435, 0, 3, 6,
                                                                       5335, 1645, 5395, 6229,
                                                                       6283, 6715, 2827, 6823,
                                                                       721, 751, 3475, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 7615, 0, 3, 6,
                                                                       5455, 1765, 5515, 6391,
                                                                       6445, 6931, 3043, 7039,
                                                                       811, 841, 3745, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 7795, 0, 3, 6,
                                                                       5515, 1795, 5575, 6445,
                                                                       6499, 7039, 3097, 7147,
                                                                       841, 871, 3835, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 7975, 0, 3, 6,
                                                                       5635, 1945, 5725, 6607,
                                                                       6715, 7255, 3385, 7435,
                                                                       931, 976, 4195, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgd_three_center_electron_repulsion_0(buffer, 8245, 0, 3, 6,
                                                                       5815, 2125, 5905, 6931,
                                                                       7039, 7615, 3745, 7795,
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

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8935, 3, 6, 8635,
                                                                       4807, 8665, 1339, 1357,
                                                                       4987, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 8995, 3, 6, 8665,
                                                                       4825, 8695, 1357, 1375,
                                                                       5023, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9055, 3, 6, 8695,
                                                                       4843, 8725, 1375, 1393,
                                                                       5059, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9115, 3, 6, 8725,
                                                                       4861, 8755, 1393, 1411,
                                                                       5095, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9175, 3, 6, 8785,
                                                                       4897, 8815, 1447, 1465,
                                                                       5131, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9235, 3, 6, 8815,
                                                                       4915, 8845, 1465, 1483,
                                                                       5167, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9295, 3, 6, 8845,
                                                                       4933, 8875, 1483, 1501,
                                                                       5203, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 9355, 3, 6, 8875,
                                                                       4951, 8905, 1501, 1519,
                                                                       5239, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9415, 3, 6, 8935,
                                                                       4987, 8995, 1555, 1585,
                                                                       5275, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9515, 3, 6, 8995,
                                                                       5023, 9055, 1585, 1615,
                                                                       5335, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9615, 3, 6, 9055,
                                                                       5059, 9115, 1615, 1645,
                                                                       5395, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9715, 3, 6, 9175,
                                                                       5131, 9235, 1705, 1735,
                                                                       5455, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9815, 3, 6, 9235,
                                                                       5167, 9295, 1735, 1765,
                                                                       5515, ncols, gamma, p,
                                                                       q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 9915, 3, 6, 9295,
                                                                       5203, 9355, 1765, 1795,
                                                                       5575, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10015, 3, 6, 9415,
                                                                       5275, 9515, 1855, 1900,
                                                                       5635, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10165, 3, 6, 9515,
                                                                       5335, 9615, 1900, 1945,
                                                                       5725, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10315, 3, 6, 9715,
                                                                       5455, 9815, 2035, 2080,
                                                                       5815, ncols, gamma, p,
                                                                       q);

                    compute_prim_sgf_three_center_electron_repulsion_0(buffer, 10465, 3, 6, 9815,
                                                                       5515, 9915, 2080, 2125,
                                                                       5905, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10615, 0, 6, 8515,
                                                                       4735, 8525, 2215, 2224,
                                                                       5995, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10645, 0, 6, 8525,
                                                                       4741, 8535, 2224, 2233,
                                                                       6013, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10675, 0, 6, 8535,
                                                                       4747, 8545, 2233, 2242,
                                                                       6031, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10705, 0, 6, 8545,
                                                                       4753, 8555, 2242, 2251,
                                                                       6049, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10735, 0, 6, 8555,
                                                                       4759, 8565, 2251, 2260,
                                                                       6067, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10765, 0, 6, 8575,
                                                                       4771, 8585, 2278, 2287,
                                                                       6085, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10795, 0, 6, 8585,
                                                                       4777, 8595, 2287, 2296,
                                                                       6103, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10825, 0, 6, 8595,
                                                                       4783, 8605, 2296, 2305,
                                                                       6121, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10855, 0, 6, 8605,
                                                                       4789, 8615, 2305, 2314,
                                                                       6139, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10885, 0, 6, 8615,
                                                                       4795, 8625, 2314, 2323,
                                                                       6157, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10915, 0, 3, 6,
                                                                       8635, 4807, 8665, 10615,
                                                                       5995, 10645, 2341, 2368,
                                                                       6175, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11005, 0, 3, 6,
                                                                       8665, 4825, 8695, 10645,
                                                                       6013, 10675, 2368, 2395,
                                                                       6229, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11095, 0, 3, 6,
                                                                       8695, 4843, 8725, 10675,
                                                                       6031, 10705, 2395, 2422,
                                                                       6283, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11185, 0, 3, 6,
                                                                       8725, 4861, 8755, 10705,
                                                                       6049, 10735, 2422, 2449,
                                                                       6337, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11275, 0, 3, 6,
                                                                       8785, 4897, 8815, 10765,
                                                                       6085, 10795, 2503, 2530,
                                                                       6391, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11365, 0, 3, 6,
                                                                       8815, 4915, 8845, 10795,
                                                                       6103, 10825, 2530, 2557,
                                                                       6445, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11455, 0, 3, 6,
                                                                       8845, 4933, 8875, 10825,
                                                                       6121, 10855, 2557, 2584,
                                                                       6499, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11545, 0, 3, 6,
                                                                       8875, 4951, 8905, 10855,
                                                                       6139, 10885, 2584, 2611,
                                                                       6553, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 11635, 0, 3, 6,
                                                                       8935, 4987, 8995, 10915,
                                                                       6175, 11005, 2665, 2719,
                                                                       6607, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 11815, 0, 3, 6,
                                                                       8995, 5023, 9055, 11005,
                                                                       6229, 11095, 2719, 2773,
                                                                       6715, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 11995, 0, 3, 6,
                                                                       9055, 5059, 9115, 11095,
                                                                       6283, 11185, 2773, 2827,
                                                                       6823, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 12175, 0, 3, 6,
                                                                       9175, 5131, 9235, 11275,
                                                                       6391, 11365, 2935, 2989,
                                                                       6931, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 12355, 0, 3, 6,
                                                                       9235, 5167, 9295, 11365,
                                                                       6445, 11455, 2989, 3043,
                                                                       7039, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 12535, 0, 3, 6,
                                                                       9295, 5203, 9355, 11455,
                                                                       6499, 11545, 3043, 3097,
                                                                       7147, ncols, gamma, p,
                                                                       q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 12715, 0, 3, 6,
                                                                       9415, 5275, 9515, 10915,
                                                                       11005, 11635, 6607, 11815,
                                                                       3205, 3295, 7255, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 13015, 0, 3, 6,
                                                                       9515, 5335, 9615, 11005,
                                                                       11095, 11815, 6715, 11995,
                                                                       3295, 3385, 7435, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 13315, 0, 3, 6,
                                                                       9715, 5455, 9815, 11275,
                                                                       11365, 12175, 6931, 12355,
                                                                       3565, 3655, 7615, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 13615, 0, 3, 6,
                                                                       9815, 5515, 9915, 11365,
                                                                       11455, 12355, 7039, 12535,
                                                                       3655, 3745, 7795, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 13915, 0, 3, 6,
                                                                       10015, 5635, 10165, 11635,
                                                                       11815, 12715, 7255, 13015,
                                                                       3925, 4060, 7975, ncols,
                                                                       gamma, p, q);

                    compute_prim_pgf_three_center_electron_repulsion_0(buffer, 14365, 0, 3, 6,
                                                                       10315, 5815, 10465, 12175,
                                                                       12355, 13315, 7615, 13615,
                                                                       4330, 4465, 8245, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 14815, 14365, 1, 150, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 14965, 14365, 1, 150, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 15115, 14365, 1, 150, ncols, alpha);

                    simdgeo::geom_s_x(buffer, 15265, 13915, 1, 150, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 15415, 13915, 1, 150, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 15565, 13915, 1, 150, ncols, alpha);

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
