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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecSFI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecPDI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecPSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSDS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSFS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecSSP.hpp"
#include "SimdTransformF.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_sfi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_sfi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 45089, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 546 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 45089, 43279, 1680, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 9, 6, 10,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 21, 6, 10,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 33, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 36, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 39, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 42, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 45, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 48, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 51, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 54, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 57, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 60, 3, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 63, 3, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 66, 3, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 69, 3, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 72, 3, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 75, 3, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 78, 3, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 81, 3, 6, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 84, 3, 6, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 87, 3, 6, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 90, 3, 6, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 93, 3, 6, 10, 11,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 99, 3, 6, 11, 12,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 105, 3, 6, 12, 13,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 111, 3, 6, 13, 14,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 117, 3, 6, 14, 15,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 123, 3, 6, 15, 16,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 129, 3, 6, 16, 17,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 135, 3, 6, 17, 18,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 141, 3, 6, 18, 19,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 147, 3, 6, 22, 23,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 153, 3, 6, 23, 24,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 159, 3, 6, 24, 25,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 165, 3, 6, 25, 26,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 171, 3, 6, 26, 27,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 177, 3, 6, 27, 28,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 183, 3, 6, 28, 29,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 189, 3, 6, 29, 30,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_sds_three_center_electron_repulsion_0(buffer, 195, 3, 6, 30, 31,
                                                                       87, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 201, 3, 6, 33, 36,
                                                                       93, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 211, 3, 6, 36, 39,
                                                                       99, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 221, 3, 6, 39, 42,
                                                                       105, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 231, 3, 6, 42, 45,
                                                                       111, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 241, 3, 6, 45, 48,
                                                                       117, 123, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 251, 3, 6, 48, 51,
                                                                       123, 129, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 261, 3, 6, 51, 54,
                                                                       129, 135, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 271, 3, 6, 54, 57,
                                                                       135, 141, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 281, 3, 6, 63, 66,
                                                                       147, 153, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 291, 3, 6, 66, 69,
                                                                       153, 159, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 301, 3, 6, 69, 72,
                                                                       159, 165, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 311, 3, 6, 72, 75,
                                                                       165, 171, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 321, 3, 6, 75, 78,
                                                                       171, 177, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 331, 3, 6, 78, 81,
                                                                       177, 183, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 341, 3, 6, 81, 84,
                                                                       183, 189, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfs_three_center_electron_repulsion_0(buffer, 351, 3, 6, 84, 87,
                                                                       189, 195, ncols, gamma, p,
                                                                       q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 361, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 364, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 367, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 370, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 373, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 376, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 379, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 382, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 385, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 388, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 391, 0, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 394, 0, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 397, 0, 6, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 400, 0, 6, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 403, 0, 6, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 406, 0, 6, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 409, 0, 6, 10, 11,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 418, 0, 6, 11, 12,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 427, 0, 6, 12, 13,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 436, 0, 6, 13, 14,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 445, 0, 6, 14, 15,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 454, 0, 6, 15, 16,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 463, 0, 6, 16, 17,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 472, 0, 6, 17, 18,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 481, 0, 6, 18, 19,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 490, 0, 6, 22, 23,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 499, 0, 6, 23, 24,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 508, 0, 6, 24, 25,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 517, 0, 6, 25, 26,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 526, 0, 6, 26, 27,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 535, 0, 6, 27, 28,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 544, 0, 6, 28, 29,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 553, 0, 6, 29, 30,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 562, 0, 6, 30, 31,
                                                                       87, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 571, 0, 3, 6, 33,
                                                                       36, 93, 99, 409, 418,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 589, 0, 3, 6, 36,
                                                                       39, 99, 105, 418, 427,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 607, 0, 3, 6, 39,
                                                                       42, 105, 111, 427, 436,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 625, 0, 3, 6, 42,
                                                                       45, 111, 117, 436, 445,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 643, 0, 3, 6, 45,
                                                                       48, 117, 123, 445, 454,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 661, 0, 3, 6, 48,
                                                                       51, 123, 129, 454, 463,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 679, 0, 3, 6, 51,
                                                                       54, 129, 135, 463, 472,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 697, 0, 3, 6, 54,
                                                                       57, 135, 141, 472, 481,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 715, 0, 3, 6, 63,
                                                                       66, 147, 153, 490, 499,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 733, 0, 3, 6, 66,
                                                                       69, 153, 159, 499, 508,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 751, 0, 3, 6, 69,
                                                                       72, 159, 165, 508, 517,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 769, 0, 3, 6, 72,
                                                                       75, 165, 171, 517, 526,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 787, 0, 3, 6, 75,
                                                                       78, 171, 177, 526, 535,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 805, 0, 3, 6, 78,
                                                                       81, 177, 183, 535, 544,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 823, 0, 3, 6, 81,
                                                                       84, 183, 189, 544, 553,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 841, 0, 3, 6, 84,
                                                                       87, 189, 195, 553, 562,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 859, 0, 3, 6, 93,
                                                                       99, 201, 211, 571, 589,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 889, 0, 3, 6, 99,
                                                                       105, 211, 221, 589, 607,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 919, 0, 3, 6, 105,
                                                                       111, 221, 231, 607, 625,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 949, 0, 3, 6, 111,
                                                                       117, 231, 241, 625, 643,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 979, 0, 3, 6, 117,
                                                                       123, 241, 251, 643, 661,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1009, 0, 3, 6,
                                                                       123, 129, 251, 261, 661,
                                                                       679, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1039, 0, 3, 6,
                                                                       129, 135, 261, 271, 679,
                                                                       697, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1069, 0, 3, 6,
                                                                       147, 153, 281, 291, 715,
                                                                       733, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1099, 0, 3, 6,
                                                                       153, 159, 291, 301, 733,
                                                                       751, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1129, 0, 3, 6,
                                                                       159, 165, 301, 311, 751,
                                                                       769, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1159, 0, 3, 6,
                                                                       165, 171, 311, 321, 769,
                                                                       787, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1189, 0, 3, 6,
                                                                       171, 177, 321, 331, 787,
                                                                       805, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1219, 0, 3, 6,
                                                                       177, 183, 331, 341, 805,
                                                                       823, ncols, gamma, p, q);

                    compute_prim_pfs_three_center_electron_repulsion_0(buffer, 1249, 0, 3, 6,
                                                                       183, 189, 341, 351, 823,
                                                                       841, ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1279, 6, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1282, 6, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1285, 6, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1288, 6, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1291, 6, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1294, 6, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1297, 6, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1300, 6, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1303, 6, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1306, 6, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1309, 6, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1312, 6, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1315, 6, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1318, 6, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1321, 6, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1324, 6, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1327, 6, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1330, 6, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1333, 6, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1342, 6, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1351, 6, 14, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1360, 6, 15, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1369, 6, 16, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1378, 6, 17, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1387, 6, 18, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1396, 6, 19, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1405, 6, 24, 69,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1414, 6, 25, 72,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1423, 6, 26, 75,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1432, 6, 27, 78,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1441, 6, 28, 81,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1450, 6, 29, 84,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1459, 6, 30, 87,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1468, 6, 31, 90,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1477, 6, 39, 105,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1495, 6, 42, 111,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1513, 6, 45, 117,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1531, 6, 48, 123,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1549, 6, 51, 129,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1567, 6, 54, 135,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1585, 6, 57, 141,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1603, 6, 69, 159,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1621, 6, 72, 165,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1639, 6, 75, 171,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1657, 6, 78, 177,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1675, 6, 81, 183,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1693, 6, 84, 189,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 1711, 6, 87, 195,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1729, 6, 105, 221,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1759, 6, 111, 231,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1789, 6, 117, 241,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1819, 6, 123, 251,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1849, 6, 129, 261,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1879, 6, 135, 271,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1909, 6, 159, 301,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1939, 6, 165, 311,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1969, 6, 171, 321,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 1999, 6, 177, 331,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2029, 6, 183, 341,
                                                                       ncols, p, q);

                    compute_prim_sfp_three_center_electron_repulsion_0(buffer, 2059, 6, 189, 351,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2089, 6, 12, 361,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2098, 6, 13, 364,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2107, 6, 14, 367,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2116, 6, 15, 370,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2125, 6, 16, 373,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2134, 6, 17, 376,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2143, 6, 18, 379,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2152, 6, 19, 382,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2161, 6, 24, 385,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2170, 6, 25, 388,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2179, 6, 26, 391,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2188, 6, 27, 394,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2197, 6, 28, 397,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2206, 6, 29, 400,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2215, 6, 30, 403,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2224, 6, 31, 406,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2233, 6, 39, 361,
                                                                       427, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2260, 6, 42, 364,
                                                                       436, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2287, 6, 45, 367,
                                                                       445, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2314, 6, 48, 370,
                                                                       454, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2341, 6, 51, 373,
                                                                       463, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2368, 6, 54, 376,
                                                                       472, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2395, 6, 57, 379,
                                                                       481, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2422, 6, 69, 385,
                                                                       508, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2449, 6, 72, 388,
                                                                       517, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2476, 6, 75, 391,
                                                                       526, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2503, 6, 78, 394,
                                                                       535, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2530, 6, 81, 397,
                                                                       544, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2557, 6, 84, 400,
                                                                       553, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2584, 6, 87, 403,
                                                                       562, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2611, 3, 6, 105,
                                                                       2233, 427, 2260, 607,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2665, 3, 6, 111,
                                                                       2260, 436, 2287, 625,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2719, 3, 6, 117,
                                                                       2287, 445, 2314, 643,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2773, 3, 6, 123,
                                                                       2314, 454, 2341, 661,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2827, 3, 6, 129,
                                                                       2341, 463, 2368, 679,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2881, 3, 6, 135,
                                                                       2368, 472, 2395, 697,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2935, 3, 6, 159,
                                                                       2422, 508, 2449, 751,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 2989, 3, 6, 165,
                                                                       2449, 517, 2476, 769,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3043, 3, 6, 171,
                                                                       2476, 526, 2503, 787,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3097, 3, 6, 177,
                                                                       2503, 535, 2530, 805,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3151, 3, 6, 183,
                                                                       2530, 544, 2557, 823,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 3205, 3, 6, 189,
                                                                       2557, 553, 2584, 841,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3259, 3, 6, 221,
                                                                       2611, 607, 2665, 919,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3349, 3, 6, 231,
                                                                       2665, 625, 2719, 949,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3439, 3, 6, 241,
                                                                       2719, 643, 2773, 979,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3529, 3, 6, 251,
                                                                       2773, 661, 2827, 1009,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3619, 3, 6, 261,
                                                                       2827, 679, 2881, 1039,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3709, 3, 6, 301,
                                                                       2935, 751, 2989, 1129,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3799, 3, 6, 311,
                                                                       2989, 769, 3043, 1159,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3889, 3, 6, 321,
                                                                       3043, 787, 3097, 1189,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 3979, 3, 6, 331,
                                                                       3097, 805, 3151, 1219,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfp_three_center_electron_repulsion_0(buffer, 4069, 3, 6, 341,
                                                                       3151, 823, 3205, 1249,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4159, 6, 10, 11,
                                                                       1279, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4165, 6, 11, 12,
                                                                       1282, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4171, 6, 12, 13,
                                                                       1285, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4177, 6, 13, 14,
                                                                       1288, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4183, 6, 14, 15,
                                                                       1291, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4189, 6, 15, 16,
                                                                       1294, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4195, 6, 16, 17,
                                                                       1297, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4201, 6, 17, 18,
                                                                       1300, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4207, 6, 18, 19,
                                                                       1303, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4213, 6, 22, 23,
                                                                       1306, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4219, 6, 23, 24,
                                                                       1309, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4225, 6, 24, 25,
                                                                       1312, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4231, 6, 25, 26,
                                                                       1315, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4237, 6, 26, 27,
                                                                       1318, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4243, 6, 27, 28,
                                                                       1321, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4249, 6, 28, 29,
                                                                       1324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4255, 6, 29, 30,
                                                                       1327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4261, 6, 30, 31,
                                                                       1330, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4267, 3, 6, 4159,
                                                                       1279, 4165, 1333, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4285, 3, 6, 4165,
                                                                       1282, 4171, 1342, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4303, 3, 6, 4171,
                                                                       1285, 4177, 1351, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4321, 3, 6, 4177,
                                                                       1288, 4183, 1360, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4339, 3, 6, 4183,
                                                                       1291, 4189, 1369, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4357, 3, 6, 4189,
                                                                       1294, 4195, 1378, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4375, 3, 6, 4195,
                                                                       1297, 4201, 1387, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4393, 3, 6, 4201,
                                                                       1300, 4207, 1396, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4411, 3, 6, 4213,
                                                                       1306, 4219, 1405, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4429, 3, 6, 4219,
                                                                       1309, 4225, 1414, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4447, 3, 6, 4225,
                                                                       1312, 4231, 1423, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4465, 3, 6, 4231,
                                                                       1315, 4237, 1432, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4483, 3, 6, 4237,
                                                                       1318, 4243, 1441, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4501, 3, 6, 4243,
                                                                       1321, 4249, 1450, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4519, 3, 6, 4249,
                                                                       1324, 4255, 1459, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4537, 3, 6, 4255,
                                                                       1327, 4261, 1468, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4555, 3, 6, 4267,
                                                                       1333, 4285, 93, 99, 1477,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4591, 3, 6, 4285,
                                                                       1342, 4303, 99, 105, 1495,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4627, 3, 6, 4303,
                                                                       1351, 4321, 105, 111,
                                                                       1513, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4663, 3, 6, 4321,
                                                                       1360, 4339, 111, 117,
                                                                       1531, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4699, 3, 6, 4339,
                                                                       1369, 4357, 117, 123,
                                                                       1549, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4735, 3, 6, 4357,
                                                                       1378, 4375, 123, 129,
                                                                       1567, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4771, 3, 6, 4375,
                                                                       1387, 4393, 129, 135,
                                                                       1585, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4807, 3, 6, 4411,
                                                                       1405, 4429, 147, 153,
                                                                       1603, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4843, 3, 6, 4429,
                                                                       1414, 4447, 153, 159,
                                                                       1621, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4879, 3, 6, 4447,
                                                                       1423, 4465, 159, 165,
                                                                       1639, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4915, 3, 6, 4465,
                                                                       1432, 4483, 165, 171,
                                                                       1657, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4951, 3, 6, 4483,
                                                                       1441, 4501, 171, 177,
                                                                       1675, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 4987, 3, 6, 4501,
                                                                       1450, 4519, 177, 183,
                                                                       1693, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 5023, 3, 6, 4519,
                                                                       1459, 4537, 183, 189,
                                                                       1711, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5059, 3, 6, 4555,
                                                                       1477, 4591, 201, 211,
                                                                       1729, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5119, 3, 6, 4591,
                                                                       1495, 4627, 211, 221,
                                                                       1759, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5179, 3, 6, 4627,
                                                                       1513, 4663, 221, 231,
                                                                       1789, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5239, 3, 6, 4663,
                                                                       1531, 4699, 231, 241,
                                                                       1819, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5299, 3, 6, 4699,
                                                                       1549, 4735, 241, 251,
                                                                       1849, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5359, 3, 6, 4735,
                                                                       1567, 4771, 251, 261,
                                                                       1879, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5419, 3, 6, 4807,
                                                                       1603, 4843, 281, 291,
                                                                       1909, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5479, 3, 6, 4843,
                                                                       1621, 4879, 291, 301,
                                                                       1939, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5539, 3, 6, 4879,
                                                                       1639, 4915, 301, 311,
                                                                       1969, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5599, 3, 6, 4915,
                                                                       1657, 4951, 311, 321,
                                                                       1999, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5659, 3, 6, 4951,
                                                                       1675, 4987, 321, 331,
                                                                       2029, ncols, gamma, p,
                                                                       q);

                    compute_prim_sfd_three_center_electron_repulsion_0(buffer, 5719, 3, 6, 4987,
                                                                       1693, 5023, 331, 341,
                                                                       2059, ncols, gamma, p,
                                                                       q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5779, 0, 6, 4159,
                                                                       1279, 4165, 2089, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5797, 0, 6, 4165,
                                                                       1282, 4171, 2098, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5815, 0, 6, 4171,
                                                                       1285, 4177, 2107, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5833, 0, 6, 4177,
                                                                       1288, 4183, 2116, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5851, 0, 6, 4183,
                                                                       1291, 4189, 2125, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5869, 0, 6, 4189,
                                                                       1294, 4195, 2134, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5887, 0, 6, 4195,
                                                                       1297, 4201, 2143, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5905, 0, 6, 4201,
                                                                       1300, 4207, 2152, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5923, 0, 6, 4213,
                                                                       1306, 4219, 2161, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5941, 0, 6, 4219,
                                                                       1309, 4225, 2170, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5959, 0, 6, 4225,
                                                                       1312, 4231, 2179, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5977, 0, 6, 4231,
                                                                       1315, 4237, 2188, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 5995, 0, 6, 4237,
                                                                       1318, 4243, 2197, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6013, 0, 6, 4243,
                                                                       1321, 4249, 2206, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6031, 0, 6, 4249,
                                                                       1324, 4255, 2215, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 6049, 0, 6, 4255,
                                                                       1327, 4261, 2224, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6067, 0, 3, 6,
                                                                       4267, 1333, 4285, 5779,
                                                                       2089, 5797, 409, 418,
                                                                       2233, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6121, 0, 3, 6,
                                                                       4285, 1342, 4303, 5797,
                                                                       2098, 5815, 418, 427,
                                                                       2260, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6175, 0, 3, 6,
                                                                       4303, 1351, 4321, 5815,
                                                                       2107, 5833, 427, 436,
                                                                       2287, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6229, 0, 3, 6,
                                                                       4321, 1360, 4339, 5833,
                                                                       2116, 5851, 436, 445,
                                                                       2314, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6283, 0, 3, 6,
                                                                       4339, 1369, 4357, 5851,
                                                                       2125, 5869, 445, 454,
                                                                       2341, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6337, 0, 3, 6,
                                                                       4357, 1378, 4375, 5869,
                                                                       2134, 5887, 454, 463,
                                                                       2368, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6391, 0, 3, 6,
                                                                       4375, 1387, 4393, 5887,
                                                                       2143, 5905, 463, 472,
                                                                       2395, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6445, 0, 3, 6,
                                                                       4411, 1405, 4429, 5923,
                                                                       2161, 5941, 490, 499,
                                                                       2422, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6499, 0, 3, 6,
                                                                       4429, 1414, 4447, 5941,
                                                                       2170, 5959, 499, 508,
                                                                       2449, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6553, 0, 3, 6,
                                                                       4447, 1423, 4465, 5959,
                                                                       2179, 5977, 508, 517,
                                                                       2476, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6607, 0, 3, 6,
                                                                       4465, 1432, 4483, 5977,
                                                                       2188, 5995, 517, 526,
                                                                       2503, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6661, 0, 3, 6,
                                                                       4483, 1441, 4501, 5995,
                                                                       2197, 6013, 526, 535,
                                                                       2530, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6715, 0, 3, 6,
                                                                       4501, 1450, 4519, 6013,
                                                                       2206, 6031, 535, 544,
                                                                       2557, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 6769, 0, 3, 6,
                                                                       4519, 1459, 4537, 6031,
                                                                       2215, 6049, 544, 553,
                                                                       2584, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6823, 0, 3, 6,
                                                                       4555, 1477, 4591, 6067,
                                                                       2233, 6121, 571, 589,
                                                                       2611, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 6931, 0, 3, 6,
                                                                       4591, 1495, 4627, 6121,
                                                                       2260, 6175, 589, 607,
                                                                       2665, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 7039, 0, 3, 6,
                                                                       4627, 1513, 4663, 6175,
                                                                       2287, 6229, 607, 625,
                                                                       2719, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 7147, 0, 3, 6,
                                                                       4663, 1531, 4699, 6229,
                                                                       2314, 6283, 625, 643,
                                                                       2773, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 7255, 0, 3, 6,
                                                                       4699, 1549, 4735, 6283,
                                                                       2341, 6337, 643, 661,
                                                                       2827, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 7363, 0, 3, 6,
                                                                       4735, 1567, 4771, 6337,
                                                                       2368, 6391, 661, 679,
                                                                       2881, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 7471, 0, 3, 6,
                                                                       4807, 1603, 4843, 6445,
                                                                       2422, 6499, 715, 733,
                                                                       2935, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 7579, 0, 3, 6,
                                                                       4843, 1621, 4879, 6499,
                                                                       2449, 6553, 733, 751,
                                                                       2989, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 7687, 0, 3, 6,
                                                                       4879, 1639, 4915, 6553,
                                                                       2476, 6607, 751, 769,
                                                                       3043, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 7795, 0, 3, 6,
                                                                       4915, 1657, 4951, 6607,
                                                                       2503, 6661, 769, 787,
                                                                       3097, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 7903, 0, 3, 6,
                                                                       4951, 1675, 4987, 6661,
                                                                       2530, 6715, 787, 805,
                                                                       3151, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 8011, 0, 3, 6,
                                                                       4987, 1693, 5023, 6715,
                                                                       2557, 6769, 805, 823,
                                                                       3205, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 8119, 0, 3, 6,
                                                                       5059, 1729, 5119, 6067,
                                                                       6121, 6823, 2611, 6931,
                                                                       859, 889, 3259, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 8299, 0, 3, 6,
                                                                       5119, 1759, 5179, 6121,
                                                                       6175, 6931, 2665, 7039,
                                                                       889, 919, 3349, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 8479, 0, 3, 6,
                                                                       5179, 1789, 5239, 6175,
                                                                       6229, 7039, 2719, 7147,
                                                                       919, 949, 3439, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 8659, 0, 3, 6,
                                                                       5239, 1819, 5299, 6229,
                                                                       6283, 7147, 2773, 7255,
                                                                       949, 979, 3529, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 8839, 0, 3, 6,
                                                                       5299, 1849, 5359, 6283,
                                                                       6337, 7255, 2827, 7363,
                                                                       979, 1009, 3619, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 9019, 0, 3, 6,
                                                                       5419, 1909, 5479, 6445,
                                                                       6499, 7471, 2935, 7579,
                                                                       1069, 1099, 3709, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 9199, 0, 3, 6,
                                                                       5479, 1939, 5539, 6499,
                                                                       6553, 7579, 2989, 7687,
                                                                       1099, 1129, 3799, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 9379, 0, 3, 6,
                                                                       5539, 1969, 5599, 6553,
                                                                       6607, 7687, 3043, 7795,
                                                                       1129, 1159, 3889, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 9559, 0, 3, 6,
                                                                       5599, 1999, 5659, 6607,
                                                                       6661, 7795, 3097, 7903,
                                                                       1159, 1189, 3979, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfd_three_center_electron_repulsion_0(buffer, 9739, 0, 3, 6,
                                                                       5659, 2029, 5719, 6661,
                                                                       6715, 7903, 3151, 8011,
                                                                       1189, 1219, 4069, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9919, 6, 1279,
                                                                       1282, 4171, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9929, 6, 1282,
                                                                       1285, 4177, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9939, 6, 1285,
                                                                       1288, 4183, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9949, 6, 1288,
                                                                       1291, 4189, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9959, 6, 1291,
                                                                       1294, 4195, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9969, 6, 1294,
                                                                       1297, 4201, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9979, 6, 1297,
                                                                       1300, 4207, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9989, 6, 1306,
                                                                       1309, 4225, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9999, 6, 1309,
                                                                       1312, 4231, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10009, 6, 1312,
                                                                       1315, 4237, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10019, 6, 1315,
                                                                       1318, 4243, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10029, 6, 1318,
                                                                       1321, 4249, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10039, 6, 1321,
                                                                       1324, 4255, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10049, 6, 1324,
                                                                       1327, 4261, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10059, 3, 6, 9919,
                                                                       4171, 9929, 4303, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10089, 3, 6, 9929,
                                                                       4177, 9939, 4321, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10119, 3, 6, 9939,
                                                                       4183, 9949, 4339, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10149, 3, 6, 9949,
                                                                       4189, 9959, 4357, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10179, 3, 6, 9959,
                                                                       4195, 9969, 4375, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10209, 3, 6, 9969,
                                                                       4201, 9979, 4393, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10239, 3, 6, 9989,
                                                                       4225, 9999, 4447, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10269, 3, 6, 9999,
                                                                       4231, 10009, 4465, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10299, 3, 6,
                                                                       10009, 4237, 10019, 4483,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10329, 3, 6,
                                                                       10019, 4243, 10029, 4501,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10359, 3, 6,
                                                                       10029, 4249, 10039, 4519,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10389, 3, 6,
                                                                       10039, 4255, 10049, 4537,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10419, 3, 6,
                                                                       10059, 4303, 10089, 1477,
                                                                       1495, 4627, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10479, 3, 6,
                                                                       10089, 4321, 10119, 1495,
                                                                       1513, 4663, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10539, 3, 6,
                                                                       10119, 4339, 10149, 1513,
                                                                       1531, 4699, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10599, 3, 6,
                                                                       10149, 4357, 10179, 1531,
                                                                       1549, 4735, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10659, 3, 6,
                                                                       10179, 4375, 10209, 1549,
                                                                       1567, 4771, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10719, 3, 6,
                                                                       10239, 4447, 10269, 1603,
                                                                       1621, 4879, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10779, 3, 6,
                                                                       10269, 4465, 10299, 1621,
                                                                       1639, 4915, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10839, 3, 6,
                                                                       10299, 4483, 10329, 1639,
                                                                       1657, 4951, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10899, 3, 6,
                                                                       10329, 4501, 10359, 1657,
                                                                       1675, 4987, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 10959, 3, 6,
                                                                       10359, 4519, 10389, 1675,
                                                                       1693, 5023, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11019, 3, 6,
                                                                       10419, 4627, 10479, 1729,
                                                                       1759, 5179, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11119, 3, 6,
                                                                       10479, 4663, 10539, 1759,
                                                                       1789, 5239, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11219, 3, 6,
                                                                       10539, 4699, 10599, 1789,
                                                                       1819, 5299, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11319, 3, 6,
                                                                       10599, 4735, 10659, 1819,
                                                                       1849, 5359, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11419, 3, 6,
                                                                       10719, 4879, 10779, 1909,
                                                                       1939, 5539, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11519, 3, 6,
                                                                       10779, 4915, 10839, 1939,
                                                                       1969, 5599, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11619, 3, 6,
                                                                       10839, 4951, 10899, 1969,
                                                                       1999, 5659, ncols, gamma,
                                                                       p, q);

                    compute_prim_sff_three_center_electron_repulsion_0(buffer, 11719, 3, 6,
                                                                       10899, 4987, 10959, 1999,
                                                                       2029, 5719, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 11819, 0, 6, 9919,
                                                                       4171, 9929, 2089, 2098,
                                                                       5815, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 11849, 0, 6, 9929,
                                                                       4177, 9939, 2098, 2107,
                                                                       5833, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 11879, 0, 6, 9939,
                                                                       4183, 9949, 2107, 2116,
                                                                       5851, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 11909, 0, 6, 9949,
                                                                       4189, 9959, 2116, 2125,
                                                                       5869, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 11939, 0, 6, 9959,
                                                                       4195, 9969, 2125, 2134,
                                                                       5887, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 11969, 0, 6, 9969,
                                                                       4201, 9979, 2134, 2143,
                                                                       5905, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 11999, 0, 6, 9989,
                                                                       4225, 9999, 2161, 2170,
                                                                       5959, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12029, 0, 6, 9999,
                                                                       4231, 10009, 2170, 2179,
                                                                       5977, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12059, 0, 6,
                                                                       10009, 4237, 10019, 2179,
                                                                       2188, 5995, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12089, 0, 6,
                                                                       10019, 4243, 10029, 2188,
                                                                       2197, 6013, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12119, 0, 6,
                                                                       10029, 4249, 10039, 2197,
                                                                       2206, 6031, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 12149, 0, 6,
                                                                       10039, 4255, 10049, 2206,
                                                                       2215, 6049, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 12179, 0, 3, 6,
                                                                       10059, 4303, 10089, 11819,
                                                                       5815, 11849, 2233, 2260,
                                                                       6175, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 12269, 0, 3, 6,
                                                                       10089, 4321, 10119, 11849,
                                                                       5833, 11879, 2260, 2287,
                                                                       6229, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 12359, 0, 3, 6,
                                                                       10119, 4339, 10149, 11879,
                                                                       5851, 11909, 2287, 2314,
                                                                       6283, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 12449, 0, 3, 6,
                                                                       10149, 4357, 10179, 11909,
                                                                       5869, 11939, 2314, 2341,
                                                                       6337, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 12539, 0, 3, 6,
                                                                       10179, 4375, 10209, 11939,
                                                                       5887, 11969, 2341, 2368,
                                                                       6391, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 12629, 0, 3, 6,
                                                                       10239, 4447, 10269, 11999,
                                                                       5959, 12029, 2422, 2449,
                                                                       6553, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 12719, 0, 3, 6,
                                                                       10269, 4465, 10299, 12029,
                                                                       5977, 12059, 2449, 2476,
                                                                       6607, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 12809, 0, 3, 6,
                                                                       10299, 4483, 10329, 12059,
                                                                       5995, 12089, 2476, 2503,
                                                                       6661, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 12899, 0, 3, 6,
                                                                       10329, 4501, 10359, 12089,
                                                                       6013, 12119, 2503, 2530,
                                                                       6715, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 12989, 0, 3, 6,
                                                                       10359, 4519, 10389, 12119,
                                                                       6031, 12149, 2530, 2557,
                                                                       6769, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 13079, 0, 3, 6,
                                                                       10419, 4627, 10479, 12179,
                                                                       6175, 12269, 2611, 2665,
                                                                       7039, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 13259, 0, 3, 6,
                                                                       10479, 4663, 10539, 12269,
                                                                       6229, 12359, 2665, 2719,
                                                                       7147, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 13439, 0, 3, 6,
                                                                       10539, 4699, 10599, 12359,
                                                                       6283, 12449, 2719, 2773,
                                                                       7255, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 13619, 0, 3, 6,
                                                                       10599, 4735, 10659, 12449,
                                                                       6337, 12539, 2773, 2827,
                                                                       7363, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 13799, 0, 3, 6,
                                                                       10719, 4879, 10779, 12629,
                                                                       6553, 12719, 2935, 2989,
                                                                       7687, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 13979, 0, 3, 6,
                                                                       10779, 4915, 10839, 12719,
                                                                       6607, 12809, 2989, 3043,
                                                                       7795, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 14159, 0, 3, 6,
                                                                       10839, 4951, 10899, 12809,
                                                                       6661, 12899, 3043, 3097,
                                                                       7903, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 14339, 0, 3, 6,
                                                                       10899, 4987, 10959, 12899,
                                                                       6715, 12989, 3097, 3151,
                                                                       8011, ncols, gamma, p,
                                                                       q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 14519, 0, 3, 6,
                                                                       11019, 5179, 11119, 12179,
                                                                       12269, 13079, 7039, 13259,
                                                                       3259, 3349, 8479, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 14819, 0, 3, 6,
                                                                       11119, 5239, 11219, 12269,
                                                                       12359, 13259, 7147, 13439,
                                                                       3349, 3439, 8659, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 15119, 0, 3, 6,
                                                                       11219, 5299, 11319, 12359,
                                                                       12449, 13439, 7255, 13619,
                                                                       3439, 3529, 8839, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 15419, 0, 3, 6,
                                                                       11419, 5539, 11519, 12629,
                                                                       12719, 13799, 7687, 13979,
                                                                       3709, 3799, 9379, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 15719, 0, 3, 6,
                                                                       11519, 5599, 11619, 12719,
                                                                       12809, 13979, 7795, 14159,
                                                                       3799, 3889, 9559, ncols,
                                                                       gamma, p, q);

                    compute_prim_pff_three_center_electron_repulsion_0(buffer, 16019, 0, 3, 6,
                                                                       11619, 5659, 11719, 12809,
                                                                       12899, 14159, 7903, 14339,
                                                                       3889, 3979, 9739, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16319, 6, 4159,
                                                                       4165, 9919, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16334, 6, 4165,
                                                                       4171, 9929, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16349, 6, 4171,
                                                                       4177, 9939, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16364, 6, 4177,
                                                                       4183, 9949, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16379, 6, 4183,
                                                                       4189, 9959, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16394, 6, 4189,
                                                                       4195, 9969, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16409, 6, 4195,
                                                                       4201, 9979, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16424, 6, 4213,
                                                                       4219, 9989, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16439, 6, 4219,
                                                                       4225, 9999, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16454, 6, 4225,
                                                                       4231, 10009, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16469, 6, 4231,
                                                                       4237, 10019, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16484, 6, 4237,
                                                                       4243, 10029, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16499, 6, 4243,
                                                                       4249, 10039, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16514, 6, 4249,
                                                                       4255, 10049, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16529, 3, 6,
                                                                       16319, 9919, 16334, 4267,
                                                                       4285, 10059, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16574, 3, 6,
                                                                       16334, 9929, 16349, 4285,
                                                                       4303, 10089, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16619, 3, 6,
                                                                       16349, 9939, 16364, 4303,
                                                                       4321, 10119, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16664, 3, 6,
                                                                       16364, 9949, 16379, 4321,
                                                                       4339, 10149, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16709, 3, 6,
                                                                       16379, 9959, 16394, 4339,
                                                                       4357, 10179, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16754, 3, 6,
                                                                       16394, 9969, 16409, 4357,
                                                                       4375, 10209, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16799, 3, 6,
                                                                       16424, 9989, 16439, 4411,
                                                                       4429, 10239, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16844, 3, 6,
                                                                       16439, 9999, 16454, 4429,
                                                                       4447, 10269, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16889, 3, 6,
                                                                       16454, 10009, 16469, 4447,
                                                                       4465, 10299, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16934, 3, 6,
                                                                       16469, 10019, 16484, 4465,
                                                                       4483, 10329, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16979, 3, 6,
                                                                       16484, 10029, 16499, 4483,
                                                                       4501, 10359, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17024, 3, 6,
                                                                       16499, 10039, 16514, 4501,
                                                                       4519, 10389, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17069, 3, 6,
                                                                       16529, 10059, 16574, 4555,
                                                                       4591, 10419, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17159, 3, 6,
                                                                       16574, 10089, 16619, 4591,
                                                                       4627, 10479, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17249, 3, 6,
                                                                       16619, 10119, 16664, 4627,
                                                                       4663, 10539, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17339, 3, 6,
                                                                       16664, 10149, 16709, 4663,
                                                                       4699, 10599, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17429, 3, 6,
                                                                       16709, 10179, 16754, 4699,
                                                                       4735, 10659, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17519, 3, 6,
                                                                       16799, 10239, 16844, 4807,
                                                                       4843, 10719, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17609, 3, 6,
                                                                       16844, 10269, 16889, 4843,
                                                                       4879, 10779, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17699, 3, 6,
                                                                       16889, 10299, 16934, 4879,
                                                                       4915, 10839, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17789, 3, 6,
                                                                       16934, 10329, 16979, 4915,
                                                                       4951, 10899, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 17879, 3, 6,
                                                                       16979, 10359, 17024, 4951,
                                                                       4987, 10959, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 17969, 3, 6,
                                                                       17069, 10419, 17159, 5059,
                                                                       5119, 11019, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 18119, 3, 6,
                                                                       17159, 10479, 17249, 5119,
                                                                       5179, 11119, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 18269, 3, 6,
                                                                       17249, 10539, 17339, 5179,
                                                                       5239, 11219, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 18419, 3, 6,
                                                                       17339, 10599, 17429, 5239,
                                                                       5299, 11319, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 18569, 3, 6,
                                                                       17519, 10719, 17609, 5419,
                                                                       5479, 11419, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 18719, 3, 6,
                                                                       17609, 10779, 17699, 5479,
                                                                       5539, 11519, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 18869, 3, 6,
                                                                       17699, 10839, 17789, 5539,
                                                                       5599, 11619, ncols, gamma,
                                                                       p, q);

                    compute_prim_sfg_three_center_electron_repulsion_0(buffer, 19019, 3, 6,
                                                                       17789, 10899, 17879, 5599,
                                                                       5659, 11719, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19169, 0, 6,
                                                                       16319, 9919, 16334, 5779,
                                                                       5797, 11819, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19214, 0, 6,
                                                                       16334, 9929, 16349, 5797,
                                                                       5815, 11849, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19259, 0, 6,
                                                                       16349, 9939, 16364, 5815,
                                                                       5833, 11879, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19304, 0, 6,
                                                                       16364, 9949, 16379, 5833,
                                                                       5851, 11909, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19349, 0, 6,
                                                                       16379, 9959, 16394, 5851,
                                                                       5869, 11939, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19394, 0, 6,
                                                                       16394, 9969, 16409, 5869,
                                                                       5887, 11969, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19439, 0, 6,
                                                                       16424, 9989, 16439, 5923,
                                                                       5941, 11999, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19484, 0, 6,
                                                                       16439, 9999, 16454, 5941,
                                                                       5959, 12029, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19529, 0, 6,
                                                                       16454, 10009, 16469, 5959,
                                                                       5977, 12059, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19574, 0, 6,
                                                                       16469, 10019, 16484, 5977,
                                                                       5995, 12089, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19619, 0, 6,
                                                                       16484, 10029, 16499, 5995,
                                                                       6013, 12119, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 19664, 0, 6,
                                                                       16499, 10039, 16514, 6013,
                                                                       6031, 12149, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 19709, 0, 3, 6,
                                                                       16529, 10059, 16574,
                                                                       19169, 11819, 19214, 6067,
                                                                       6121, 12179, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 19844, 0, 3, 6,
                                                                       16574, 10089, 16619,
                                                                       19214, 11849, 19259, 6121,
                                                                       6175, 12269, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 19979, 0, 3, 6,
                                                                       16619, 10119, 16664,
                                                                       19259, 11879, 19304, 6175,
                                                                       6229, 12359, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 20114, 0, 3, 6,
                                                                       16664, 10149, 16709,
                                                                       19304, 11909, 19349, 6229,
                                                                       6283, 12449, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 20249, 0, 3, 6,
                                                                       16709, 10179, 16754,
                                                                       19349, 11939, 19394, 6283,
                                                                       6337, 12539, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 20384, 0, 3, 6,
                                                                       16799, 10239, 16844,
                                                                       19439, 11999, 19484, 6445,
                                                                       6499, 12629, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 20519, 0, 3, 6,
                                                                       16844, 10269, 16889,
                                                                       19484, 12029, 19529, 6499,
                                                                       6553, 12719, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 20654, 0, 3, 6,
                                                                       16889, 10299, 16934,
                                                                       19529, 12059, 19574, 6553,
                                                                       6607, 12809, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 20789, 0, 3, 6,
                                                                       16934, 10329, 16979,
                                                                       19574, 12089, 19619, 6607,
                                                                       6661, 12899, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 20924, 0, 3, 6,
                                                                       16979, 10359, 17024,
                                                                       19619, 12119, 19664, 6661,
                                                                       6715, 12989, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 21059, 0, 3, 6,
                                                                       17069, 10419, 17159,
                                                                       19709, 12179, 19844, 6823,
                                                                       6931, 13079, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 21329, 0, 3, 6,
                                                                       17159, 10479, 17249,
                                                                       19844, 12269, 19979, 6931,
                                                                       7039, 13259, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 21599, 0, 3, 6,
                                                                       17249, 10539, 17339,
                                                                       19979, 12359, 20114, 7039,
                                                                       7147, 13439, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 21869, 0, 3, 6,
                                                                       17339, 10599, 17429,
                                                                       20114, 12449, 20249, 7147,
                                                                       7255, 13619, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 22139, 0, 3, 6,
                                                                       17519, 10719, 17609,
                                                                       20384, 12629, 20519, 7471,
                                                                       7579, 13799, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 22409, 0, 3, 6,
                                                                       17609, 10779, 17699,
                                                                       20519, 12719, 20654, 7579,
                                                                       7687, 13979, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 22679, 0, 3, 6,
                                                                       17699, 10839, 17789,
                                                                       20654, 12809, 20789, 7687,
                                                                       7795, 14159, ncols, gamma,
                                                                       p, q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 22949, 0, 3, 6,
                                                                       17789, 10899, 17879,
                                                                       20789, 12899, 20924, 7795,
                                                                       7903, 14339, ncols, gamma,
                                                                       p, q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 23219, 0, 3, 6,
                                                                       17969, 11019, 18119,
                                                                       19709, 19844, 21059,
                                                                       13079, 21329, 8119, 8299,
                                                                       14519, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 23669, 0, 3, 6,
                                                                       18119, 11119, 18269,
                                                                       19844, 19979, 21329,
                                                                       13259, 21599, 8299, 8479,
                                                                       14819, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 24119, 0, 3, 6,
                                                                       18269, 11219, 18419,
                                                                       19979, 20114, 21599,
                                                                       13439, 21869, 8479, 8659,
                                                                       15119, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 24569, 0, 3, 6,
                                                                       18569, 11419, 18719,
                                                                       20384, 20519, 22139,
                                                                       13799, 22409, 9019, 9199,
                                                                       15419, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 25019, 0, 3, 6,
                                                                       18719, 11519, 18869,
                                                                       20519, 20654, 22409,
                                                                       13979, 22679, 9199, 9379,
                                                                       15719, ncols, gamma, p,
                                                                       q);

                    compute_prim_pfg_three_center_electron_repulsion_0(buffer, 25469, 0, 3, 6,
                                                                       18869, 11619, 19019,
                                                                       20654, 20789, 22679,
                                                                       14159, 22949, 9379, 9559,
                                                                       16019, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25919, 6, 9919,
                                                                       9929, 16349, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25940, 6, 9929,
                                                                       9939, 16364, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25961, 6, 9939,
                                                                       9949, 16379, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25982, 6, 9949,
                                                                       9959, 16394, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26003, 6, 9959,
                                                                       9969, 16409, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26024, 6, 9989,
                                                                       9999, 16454, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26045, 6, 9999,
                                                                       10009, 16469, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26066, 6, 10009,
                                                                       10019, 16484, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26087, 6, 10019,
                                                                       10029, 16499, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26108, 6, 10029,
                                                                       10039, 16514, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26129, 3, 6,
                                                                       25919, 16349, 25940,
                                                                       10059, 10089, 16619,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26192, 3, 6,
                                                                       25940, 16364, 25961,
                                                                       10089, 10119, 16664,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26255, 3, 6,
                                                                       25961, 16379, 25982,
                                                                       10119, 10149, 16709,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26318, 3, 6,
                                                                       25982, 16394, 26003,
                                                                       10149, 10179, 16754,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26381, 3, 6,
                                                                       26024, 16454, 26045,
                                                                       10239, 10269, 16889,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26444, 3, 6,
                                                                       26045, 16469, 26066,
                                                                       10269, 10299, 16934,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26507, 3, 6,
                                                                       26066, 16484, 26087,
                                                                       10299, 10329, 16979,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26570, 3, 6,
                                                                       26087, 16499, 26108,
                                                                       10329, 10359, 17024,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26633, 3, 6,
                                                                       26129, 16619, 26192,
                                                                       10419, 10479, 17249,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26759, 3, 6,
                                                                       26192, 16664, 26255,
                                                                       10479, 10539, 17339,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 26885, 3, 6,
                                                                       26255, 16709, 26318,
                                                                       10539, 10599, 17429,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 27011, 3, 6,
                                                                       26381, 16889, 26444,
                                                                       10719, 10779, 17699,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 27137, 3, 6,
                                                                       26444, 16934, 26507,
                                                                       10779, 10839, 17789,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 27263, 3, 6,
                                                                       26507, 16979, 26570,
                                                                       10839, 10899, 17879,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 27389, 3, 6,
                                                                       26633, 17249, 26759,
                                                                       11019, 11119, 18269,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 27599, 3, 6,
                                                                       26759, 17339, 26885,
                                                                       11119, 11219, 18419,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 27809, 3, 6,
                                                                       27011, 17699, 27137,
                                                                       11419, 11519, 18869,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfh_three_center_electron_repulsion_0(buffer, 28019, 3, 6,
                                                                       27137, 17789, 27263,
                                                                       11519, 11619, 19019,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 28229, 0, 6,
                                                                       25919, 16349, 25940,
                                                                       11819, 11849, 19259,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 28292, 0, 6,
                                                                       25940, 16364, 25961,
                                                                       11849, 11879, 19304,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 28355, 0, 6,
                                                                       25961, 16379, 25982,
                                                                       11879, 11909, 19349,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 28418, 0, 6,
                                                                       25982, 16394, 26003,
                                                                       11909, 11939, 19394,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 28481, 0, 6,
                                                                       26024, 16454, 26045,
                                                                       11999, 12029, 19529,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 28544, 0, 6,
                                                                       26045, 16469, 26066,
                                                                       12029, 12059, 19574,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 28607, 0, 6,
                                                                       26066, 16484, 26087,
                                                                       12059, 12089, 19619,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 28670, 0, 6,
                                                                       26087, 16499, 26108,
                                                                       12089, 12119, 19664,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 28733, 0, 3, 6,
                                                                       26129, 16619, 26192,
                                                                       28229, 19259, 28292,
                                                                       12179, 12269, 19979,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 28922, 0, 3, 6,
                                                                       26192, 16664, 26255,
                                                                       28292, 19304, 28355,
                                                                       12269, 12359, 20114,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 29111, 0, 3, 6,
                                                                       26255, 16709, 26318,
                                                                       28355, 19349, 28418,
                                                                       12359, 12449, 20249,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 29300, 0, 3, 6,
                                                                       26381, 16889, 26444,
                                                                       28481, 19529, 28544,
                                                                       12629, 12719, 20654,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 29489, 0, 3, 6,
                                                                       26444, 16934, 26507,
                                                                       28544, 19574, 28607,
                                                                       12719, 12809, 20789,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 29678, 0, 3, 6,
                                                                       26507, 16979, 26570,
                                                                       28607, 19619, 28670,
                                                                       12809, 12899, 20924,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 29867, 0, 3, 6,
                                                                       26633, 17249, 26759,
                                                                       28733, 19979, 28922,
                                                                       13079, 13259, 21599,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 30245, 0, 3, 6,
                                                                       26759, 17339, 26885,
                                                                       28922, 20114, 29111,
                                                                       13259, 13439, 21869,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 30623, 0, 3, 6,
                                                                       27011, 17699, 27137,
                                                                       29300, 20654, 29489,
                                                                       13799, 13979, 22679,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 31001, 0, 3, 6,
                                                                       27137, 17789, 27263,
                                                                       29489, 20789, 29678,
                                                                       13979, 14159, 22949,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 31379, 0, 3, 6,
                                                                       27389, 18269, 27599,
                                                                       28733, 28922, 29867,
                                                                       21599, 30245, 14519,
                                                                       14819, 24119, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfh_three_center_electron_repulsion_0(buffer, 32009, 0, 3, 6,
                                                                       27809, 18869, 28019,
                                                                       29300, 29489, 30623,
                                                                       22679, 31001, 15419,
                                                                       15719, 25469, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32639, 6, 16319,
                                                                       16334, 25919, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32667, 6, 16334,
                                                                       16349, 25940, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32695, 6, 16349,
                                                                       16364, 25961, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32723, 6, 16364,
                                                                       16379, 25982, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32751, 6, 16379,
                                                                       16394, 26003, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32779, 6, 16424,
                                                                       16439, 26024, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32807, 6, 16439,
                                                                       16454, 26045, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32835, 6, 16454,
                                                                       16469, 26066, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32863, 6, 16469,
                                                                       16484, 26087, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32891, 6, 16484,
                                                                       16499, 26108, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 32919, 3, 6,
                                                                       32639, 25919, 32667,
                                                                       16529, 16574, 26129,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33003, 3, 6,
                                                                       32667, 25940, 32695,
                                                                       16574, 16619, 26192,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33087, 3, 6,
                                                                       32695, 25961, 32723,
                                                                       16619, 16664, 26255,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33171, 3, 6,
                                                                       32723, 25982, 32751,
                                                                       16664, 16709, 26318,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33255, 3, 6,
                                                                       32779, 26024, 32807,
                                                                       16799, 16844, 26381,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33339, 3, 6,
                                                                       32807, 26045, 32835,
                                                                       16844, 16889, 26444,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33423, 3, 6,
                                                                       32835, 26066, 32863,
                                                                       16889, 16934, 26507,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33507, 3, 6,
                                                                       32863, 26087, 32891,
                                                                       16934, 16979, 26570,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 33591, 3, 6,
                                                                       32919, 26129, 33003,
                                                                       17069, 17159, 26633,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 33759, 3, 6,
                                                                       33003, 26192, 33087,
                                                                       17159, 17249, 26759,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 33927, 3, 6,
                                                                       33087, 26255, 33171,
                                                                       17249, 17339, 26885,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 34095, 3, 6,
                                                                       33255, 26381, 33339,
                                                                       17519, 17609, 27011,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 34263, 3, 6,
                                                                       33339, 26444, 33423,
                                                                       17609, 17699, 27137,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 34431, 3, 6,
                                                                       33423, 26507, 33507,
                                                                       17699, 17789, 27263,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 34599, 3, 6,
                                                                       33591, 26633, 33759,
                                                                       17969, 18119, 27389,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 34879, 3, 6,
                                                                       33759, 26759, 33927,
                                                                       18119, 18269, 27599,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 35159, 3, 6,
                                                                       34095, 27011, 34263,
                                                                       18569, 18719, 27809,
                                                                       ncols, gamma, p, q);

                    compute_prim_sfi_three_center_electron_repulsion_0(buffer, 35439, 3, 6,
                                                                       34263, 27137, 34431,
                                                                       18719, 18869, 28019,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 35719, 0, 6,
                                                                       32639, 25919, 32667,
                                                                       19169, 19214, 28229,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 35803, 0, 6,
                                                                       32667, 25940, 32695,
                                                                       19214, 19259, 28292,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 35887, 0, 6,
                                                                       32695, 25961, 32723,
                                                                       19259, 19304, 28355,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 35971, 0, 6,
                                                                       32723, 25982, 32751,
                                                                       19304, 19349, 28418,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 36055, 0, 6,
                                                                       32779, 26024, 32807,
                                                                       19439, 19484, 28481,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 36139, 0, 6,
                                                                       32807, 26045, 32835,
                                                                       19484, 19529, 28544,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 36223, 0, 6,
                                                                       32835, 26066, 32863,
                                                                       19529, 19574, 28607,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 36307, 0, 6,
                                                                       32863, 26087, 32891,
                                                                       19574, 19619, 28670,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 36391, 0, 3, 6,
                                                                       32919, 26129, 33003,
                                                                       35719, 28229, 35803,
                                                                       19709, 19844, 28733,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 36643, 0, 3, 6,
                                                                       33003, 26192, 33087,
                                                                       35803, 28292, 35887,
                                                                       19844, 19979, 28922,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 36895, 0, 3, 6,
                                                                       33087, 26255, 33171,
                                                                       35887, 28355, 35971,
                                                                       19979, 20114, 29111,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 37147, 0, 3, 6,
                                                                       33255, 26381, 33339,
                                                                       36055, 28481, 36139,
                                                                       20384, 20519, 29300,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 37399, 0, 3, 6,
                                                                       33339, 26444, 33423,
                                                                       36139, 28544, 36223,
                                                                       20519, 20654, 29489,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 37651, 0, 3, 6,
                                                                       33423, 26507, 33507,
                                                                       36223, 28607, 36307,
                                                                       20654, 20789, 29678,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 37903, 0, 3, 6,
                                                                       33591, 26633, 33759,
                                                                       36391, 28733, 36643,
                                                                       21059, 21329, 29867,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 38407, 0, 3, 6,
                                                                       33759, 26759, 33927,
                                                                       36643, 28922, 36895,
                                                                       21329, 21599, 30245,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 38911, 0, 3, 6,
                                                                       34095, 27011, 34263,
                                                                       37147, 29300, 37399,
                                                                       22139, 22409, 30623,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 39415, 0, 3, 6,
                                                                       34263, 27137, 34431,
                                                                       37399, 29489, 37651,
                                                                       22409, 22679, 31001,
                                                                       ncols, gamma, p, q);

                    compute_prim_pfi_three_center_electron_repulsion_0(buffer, 39919, 0, 3, 6,
                                                                       34599, 27389, 34879,
                                                                       36391, 36643, 37903,
                                                                       29867, 38407, 23219,
                                                                       23669, 31379, ncols,
                                                                       gamma, p, q);

                    compute_prim_pfi_three_center_electron_repulsion_0(buffer, 40759, 0, 3, 6,
                                                                       35159, 27809, 35439,
                                                                       37147, 37399, 38911,
                                                                       30623, 39415, 24569,
                                                                       25019, 32009, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 41599, 40759, 1, 280, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 41879, 40759, 1, 280, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 42159, 40759, 1, 280, ncols, alpha);

                    simdgeo::geom_s_x(buffer, 42439, 39919, 1, 280, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 42719, 39919, 1, 280, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 42999, 39919, 1, 280, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 43279, 41599, 1680, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 44959, 43279, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 44959, 13, nmax);

        simdtrf::transform_i_inner(buffer, 44959, 43559, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 91 * nvalues + n * npairs, nvalues, buffer, 44959,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 44959, 43839, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 182 * nvalues + n * npairs, nvalues, buffer, 44959,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 44959, 44119, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 273 * nvalues + n * npairs, nvalues, buffer, 44959,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 44959, 44399, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 364 * nvalues + n * npairs, nvalues, buffer, 44959,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 44959, 44679, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 455 * nvalues + n * npairs, nvalues, buffer, 44959,
                                   13, nmax);
    }

    for (size_t m = 0; m < 546; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
