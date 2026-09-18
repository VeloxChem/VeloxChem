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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecFSI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecDPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecDSS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecFSS.hpp"
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
compute_rs_geom_010_fsi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_fsi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 45101, 0, 0, dimensions);

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
        simdfunc::prepare_buffer(buffer, 45101, 43291, 1680, dimensions);

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

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 93, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 96, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 99, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 102, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 105, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 108, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 111, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 114, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 117, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 120, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 123, 0, 6, 22, 23,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 126, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 129, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 132, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 135, 0, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 138, 0, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 141, 0, 6, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 144, 0, 6, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 147, 0, 6, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 150, 0, 6, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 153, 0, 6, 10, 11,
                                                                       33, 36, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 162, 0, 6, 11, 12,
                                                                       36, 39, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 171, 0, 6, 12, 13,
                                                                       39, 42, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 180, 0, 6, 13, 14,
                                                                       42, 45, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 189, 0, 6, 14, 15,
                                                                       45, 48, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 198, 0, 6, 15, 16,
                                                                       48, 51, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 207, 0, 6, 16, 17,
                                                                       51, 54, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 216, 0, 6, 17, 18,
                                                                       54, 57, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 225, 0, 6, 18, 19,
                                                                       57, 60, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 234, 0, 6, 22, 23,
                                                                       63, 66, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 243, 0, 6, 23, 24,
                                                                       66, 69, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 252, 0, 6, 24, 25,
                                                                       69, 72, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 261, 0, 6, 25, 26,
                                                                       72, 75, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 270, 0, 6, 26, 27,
                                                                       75, 78, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 279, 0, 6, 27, 28,
                                                                       78, 81, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 288, 0, 6, 28, 29,
                                                                       81, 84, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 297, 0, 6, 29, 30,
                                                                       84, 87, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 306, 0, 6, 30, 31,
                                                                       87, 90, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 315, 0, 6, 10, 11,
                                                                       93, 96, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 321, 0, 6, 11, 12,
                                                                       96, 99, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 327, 0, 6, 12, 13,
                                                                       99, 102, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 333, 0, 6, 13, 14,
                                                                       102, 105, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 339, 0, 6, 14, 15,
                                                                       105, 108, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 345, 0, 6, 15, 16,
                                                                       108, 111, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 351, 0, 6, 16, 17,
                                                                       111, 114, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 357, 0, 6, 17, 18,
                                                                       114, 117, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 363, 0, 6, 18, 19,
                                                                       117, 120, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 369, 0, 6, 22, 23,
                                                                       123, 126, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 375, 0, 6, 23, 24,
                                                                       126, 129, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 381, 0, 6, 24, 25,
                                                                       129, 132, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 387, 0, 6, 25, 26,
                                                                       132, 135, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 393, 0, 6, 26, 27,
                                                                       135, 138, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 399, 0, 6, 27, 28,
                                                                       138, 141, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 405, 0, 6, 28, 29,
                                                                       141, 144, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 411, 0, 6, 29, 30,
                                                                       144, 147, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 417, 0, 6, 30, 31,
                                                                       147, 150, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 423, 0, 3, 6, 93,
                                                                       96, 153, 162, 315, 321,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 441, 0, 3, 6, 96,
                                                                       99, 162, 171, 321, 327,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 459, 0, 3, 6, 99,
                                                                       102, 171, 180, 327, 333,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 477, 0, 3, 6, 102,
                                                                       105, 180, 189, 333, 339,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 495, 0, 3, 6, 105,
                                                                       108, 189, 198, 339, 345,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 513, 0, 3, 6, 108,
                                                                       111, 198, 207, 345, 351,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 531, 0, 3, 6, 111,
                                                                       114, 207, 216, 351, 357,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 549, 0, 3, 6, 114,
                                                                       117, 216, 225, 357, 363,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 567, 0, 3, 6, 123,
                                                                       126, 234, 243, 369, 375,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 585, 0, 3, 6, 126,
                                                                       129, 243, 252, 375, 381,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 603, 0, 3, 6, 129,
                                                                       132, 252, 261, 381, 387,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 621, 0, 3, 6, 132,
                                                                       135, 261, 270, 387, 393,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 639, 0, 3, 6, 135,
                                                                       138, 270, 279, 393, 399,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 657, 0, 3, 6, 138,
                                                                       141, 279, 288, 399, 405,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 675, 0, 3, 6, 141,
                                                                       144, 288, 297, 405, 411,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 693, 0, 3, 6, 144,
                                                                       147, 297, 306, 411, 417,
                                                                       ncols, gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 711, 0, 6, 93, 96,
                                                                       315, 321, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 721, 0, 6, 96, 99,
                                                                       321, 327, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 731, 0, 6, 99,
                                                                       102, 327, 333, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 741, 0, 6, 102,
                                                                       105, 333, 339, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 751, 0, 6, 105,
                                                                       108, 339, 345, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 761, 0, 6, 108,
                                                                       111, 345, 351, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 771, 0, 6, 111,
                                                                       114, 351, 357, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 781, 0, 6, 114,
                                                                       117, 357, 363, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 791, 0, 6, 123,
                                                                       126, 369, 375, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 801, 0, 6, 126,
                                                                       129, 375, 381, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 811, 0, 6, 129,
                                                                       132, 381, 387, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 821, 0, 6, 132,
                                                                       135, 387, 393, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 831, 0, 6, 135,
                                                                       138, 393, 399, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 841, 0, 6, 138,
                                                                       141, 399, 405, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 851, 0, 6, 141,
                                                                       144, 405, 411, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 861, 0, 6, 144,
                                                                       147, 411, 417, ncols,
                                                                       gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 871, 0, 3, 6, 315,
                                                                       321, 423, 441, 711, 721,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 901, 0, 3, 6, 321,
                                                                       327, 441, 459, 721, 731,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 931, 0, 3, 6, 327,
                                                                       333, 459, 477, 731, 741,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 961, 0, 3, 6, 333,
                                                                       339, 477, 495, 741, 751,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 991, 0, 3, 6, 339,
                                                                       345, 495, 513, 751, 761,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1021, 0, 3, 6,
                                                                       345, 351, 513, 531, 761,
                                                                       771, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1051, 0, 3, 6,
                                                                       351, 357, 531, 549, 771,
                                                                       781, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1081, 0, 3, 6,
                                                                       369, 375, 567, 585, 791,
                                                                       801, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1111, 0, 3, 6,
                                                                       375, 381, 585, 603, 801,
                                                                       811, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1141, 0, 3, 6,
                                                                       381, 387, 603, 621, 811,
                                                                       821, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1171, 0, 3, 6,
                                                                       387, 393, 621, 639, 821,
                                                                       831, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1201, 0, 3, 6,
                                                                       393, 399, 639, 657, 831,
                                                                       841, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1231, 0, 3, 6,
                                                                       399, 405, 657, 675, 841,
                                                                       851, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1261, 0, 3, 6,
                                                                       405, 411, 675, 693, 851,
                                                                       861, ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1291, 6, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1294, 6, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1297, 6, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1300, 6, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1303, 6, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1306, 6, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1309, 6, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1312, 6, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1315, 6, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1318, 6, 24,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1321, 6, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1324, 6, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1327, 6, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1330, 6, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1333, 6, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1336, 6, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1339, 6, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1342, 6, 32,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1345, 6, 12, 39,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1354, 6, 13, 42,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1363, 6, 14, 45,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1372, 6, 15, 48,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1381, 6, 16, 51,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1390, 6, 17, 54,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1399, 6, 18, 57,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1408, 6, 19, 60,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1417, 6, 24, 69,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1426, 6, 25, 72,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1435, 6, 26, 75,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1444, 6, 27, 78,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1453, 6, 28, 81,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1462, 6, 29, 84,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1471, 6, 30, 87,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1480, 6, 31, 90,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1489, 6, 12, 99,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1498, 6, 13, 102,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1507, 6, 14, 105,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1516, 6, 15, 108,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1525, 6, 16, 111,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1534, 6, 17, 114,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1543, 6, 18, 117,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1552, 6, 19, 120,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1561, 6, 24, 129,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1570, 6, 25, 132,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1579, 6, 26, 135,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1588, 6, 27, 138,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1597, 6, 28, 141,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1606, 6, 29, 144,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1615, 6, 30, 147,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1624, 6, 31, 150,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1633, 6, 39, 99,
                                                                       171, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1660, 6, 42, 102,
                                                                       180, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1687, 6, 45, 105,
                                                                       189, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1714, 6, 48, 108,
                                                                       198, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1741, 6, 51, 111,
                                                                       207, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1768, 6, 54, 114,
                                                                       216, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1795, 6, 57, 117,
                                                                       225, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1822, 6, 69, 129,
                                                                       252, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1849, 6, 72, 132,
                                                                       261, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1876, 6, 75, 135,
                                                                       270, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1903, 6, 78, 138,
                                                                       279, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1930, 6, 81, 141,
                                                                       288, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1957, 6, 84, 144,
                                                                       297, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1984, 6, 87, 147,
                                                                       306, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2011, 6, 99, 327,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2029, 6, 102, 333,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2047, 6, 105, 339,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2065, 6, 108, 345,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2083, 6, 111, 351,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2101, 6, 114, 357,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2119, 6, 117, 363,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2137, 6, 129, 381,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2155, 6, 132, 387,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2173, 6, 135, 393,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2191, 6, 138, 399,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2209, 6, 141, 405,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2227, 6, 144, 411,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 2245, 6, 147, 417,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2263, 0, 6, 1633,
                                                                       171, 1660, 327, 459,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2317, 0, 6, 1660,
                                                                       180, 1687, 333, 477,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2371, 0, 6, 1687,
                                                                       189, 1714, 339, 495,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2425, 0, 6, 1714,
                                                                       198, 1741, 345, 513,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2479, 0, 6, 1741,
                                                                       207, 1768, 351, 531,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2533, 0, 6, 1768,
                                                                       216, 1795, 357, 549,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2587, 0, 6, 1822,
                                                                       252, 1849, 381, 603,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2641, 0, 6, 1849,
                                                                       261, 1876, 387, 621,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2695, 0, 6, 1876,
                                                                       270, 1903, 393, 639,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2749, 0, 6, 1903,
                                                                       279, 1930, 399, 657,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2803, 0, 6, 1930,
                                                                       288, 1957, 405, 675,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2857, 0, 6, 1957,
                                                                       297, 1984, 411, 693,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2911, 6, 327, 731,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2941, 6, 333, 741,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2971, 6, 339, 751,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3001, 6, 345, 761,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3031, 6, 351, 771,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3061, 6, 357, 781,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3091, 6, 381, 811,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3121, 6, 387, 821,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3151, 6, 393, 831,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3181, 6, 399, 841,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3211, 6, 405, 851,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 3241, 6, 411, 861,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3271, 0, 6, 2263,
                                                                       459, 2317, 731, 931,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3361, 0, 6, 2317,
                                                                       477, 2371, 741, 961,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3451, 0, 6, 2371,
                                                                       495, 2425, 751, 991,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3541, 0, 6, 2425,
                                                                       513, 2479, 761, 1021,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3631, 0, 6, 2479,
                                                                       531, 2533, 771, 1051,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3721, 0, 6, 2587,
                                                                       603, 2641, 811, 1141,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3811, 0, 6, 2641,
                                                                       621, 2695, 821, 1171,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3901, 0, 6, 2695,
                                                                       639, 2749, 831, 1201,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 3991, 0, 6, 2749,
                                                                       657, 2803, 841, 1231,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4081, 0, 6, 2803,
                                                                       675, 2857, 851, 1261,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4171, 6, 10, 11,
                                                                       1291, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4177, 6, 11, 12,
                                                                       1294, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4183, 6, 12, 13,
                                                                       1297, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4189, 6, 13, 14,
                                                                       1300, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4195, 6, 14, 15,
                                                                       1303, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4201, 6, 15, 16,
                                                                       1306, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4207, 6, 16, 17,
                                                                       1309, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4213, 6, 17, 18,
                                                                       1312, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4219, 6, 18, 19,
                                                                       1315, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4225, 6, 22, 23,
                                                                       1318, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4231, 6, 23, 24,
                                                                       1321, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4237, 6, 24, 25,
                                                                       1324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4243, 6, 25, 26,
                                                                       1327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4249, 6, 26, 27,
                                                                       1330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4255, 6, 27, 28,
                                                                       1333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4261, 6, 28, 29,
                                                                       1336, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4267, 6, 29, 30,
                                                                       1339, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 4273, 6, 30, 31,
                                                                       1342, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4279, 3, 6, 4171,
                                                                       1291, 4177, 1345, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4297, 3, 6, 4177,
                                                                       1294, 4183, 1354, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4315, 3, 6, 4183,
                                                                       1297, 4189, 1363, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4333, 3, 6, 4189,
                                                                       1300, 4195, 1372, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4351, 3, 6, 4195,
                                                                       1303, 4201, 1381, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4369, 3, 6, 4201,
                                                                       1306, 4207, 1390, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4387, 3, 6, 4207,
                                                                       1309, 4213, 1399, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4405, 3, 6, 4213,
                                                                       1312, 4219, 1408, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4423, 3, 6, 4225,
                                                                       1318, 4231, 1417, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4441, 3, 6, 4231,
                                                                       1321, 4237, 1426, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4459, 3, 6, 4237,
                                                                       1324, 4243, 1435, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4477, 3, 6, 4243,
                                                                       1327, 4249, 1444, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4495, 3, 6, 4249,
                                                                       1330, 4255, 1453, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4513, 3, 6, 4255,
                                                                       1333, 4261, 1462, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4531, 3, 6, 4261,
                                                                       1336, 4267, 1471, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 4549, 3, 6, 4267,
                                                                       1339, 4273, 1480, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4567, 0, 6, 4171,
                                                                       1291, 4177, 1489, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4585, 0, 6, 4177,
                                                                       1294, 4183, 1498, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4603, 0, 6, 4183,
                                                                       1297, 4189, 1507, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4621, 0, 6, 4189,
                                                                       1300, 4195, 1516, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4639, 0, 6, 4195,
                                                                       1303, 4201, 1525, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4657, 0, 6, 4201,
                                                                       1306, 4207, 1534, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4675, 0, 6, 4207,
                                                                       1309, 4213, 1543, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4693, 0, 6, 4213,
                                                                       1312, 4219, 1552, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4711, 0, 6, 4225,
                                                                       1318, 4231, 1561, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4729, 0, 6, 4231,
                                                                       1321, 4237, 1570, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4747, 0, 6, 4237,
                                                                       1324, 4243, 1579, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4765, 0, 6, 4243,
                                                                       1327, 4249, 1588, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4783, 0, 6, 4249,
                                                                       1330, 4255, 1597, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4801, 0, 6, 4255,
                                                                       1333, 4261, 1606, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4819, 0, 6, 4261,
                                                                       1336, 4267, 1615, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4837, 0, 6, 4267,
                                                                       1339, 4273, 1624, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4855, 0, 3, 6,
                                                                       4279, 1345, 4297, 4567,
                                                                       1489, 4585, 153, 162,
                                                                       1633, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4909, 0, 3, 6,
                                                                       4297, 1354, 4315, 4585,
                                                                       1498, 4603, 162, 171,
                                                                       1660, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4963, 0, 3, 6,
                                                                       4315, 1363, 4333, 4603,
                                                                       1507, 4621, 171, 180,
                                                                       1687, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5017, 0, 3, 6,
                                                                       4333, 1372, 4351, 4621,
                                                                       1516, 4639, 180, 189,
                                                                       1714, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5071, 0, 3, 6,
                                                                       4351, 1381, 4369, 4639,
                                                                       1525, 4657, 189, 198,
                                                                       1741, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5125, 0, 3, 6,
                                                                       4369, 1390, 4387, 4657,
                                                                       1534, 4675, 198, 207,
                                                                       1768, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5179, 0, 3, 6,
                                                                       4387, 1399, 4405, 4675,
                                                                       1543, 4693, 207, 216,
                                                                       1795, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5233, 0, 3, 6,
                                                                       4423, 1417, 4441, 4711,
                                                                       1561, 4729, 234, 243,
                                                                       1822, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5287, 0, 3, 6,
                                                                       4441, 1426, 4459, 4729,
                                                                       1570, 4747, 243, 252,
                                                                       1849, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5341, 0, 3, 6,
                                                                       4459, 1435, 4477, 4747,
                                                                       1579, 4765, 252, 261,
                                                                       1876, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5395, 0, 3, 6,
                                                                       4477, 1444, 4495, 4765,
                                                                       1588, 4783, 261, 270,
                                                                       1903, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5449, 0, 3, 6,
                                                                       4495, 1453, 4513, 4783,
                                                                       1597, 4801, 270, 279,
                                                                       1930, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5503, 0, 3, 6,
                                                                       4513, 1462, 4531, 4801,
                                                                       1606, 4819, 279, 288,
                                                                       1957, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 5557, 0, 3, 6,
                                                                       4531, 1471, 4549, 4819,
                                                                       1615, 4837, 288, 297,
                                                                       1984, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5611, 0, 6, 4567,
                                                                       1489, 4585, 315, 321,
                                                                       2011, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5647, 0, 6, 4585,
                                                                       1498, 4603, 321, 327,
                                                                       2029, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5683, 0, 6, 4603,
                                                                       1507, 4621, 327, 333,
                                                                       2047, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5719, 0, 6, 4621,
                                                                       1516, 4639, 333, 339,
                                                                       2065, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5755, 0, 6, 4639,
                                                                       1525, 4657, 339, 345,
                                                                       2083, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5791, 0, 6, 4657,
                                                                       1534, 4675, 345, 351,
                                                                       2101, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5827, 0, 6, 4675,
                                                                       1543, 4693, 351, 357,
                                                                       2119, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5863, 0, 6, 4711,
                                                                       1561, 4729, 369, 375,
                                                                       2137, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5899, 0, 6, 4729,
                                                                       1570, 4747, 375, 381,
                                                                       2155, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5935, 0, 6, 4747,
                                                                       1579, 4765, 381, 387,
                                                                       2173, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 5971, 0, 6, 4765,
                                                                       1588, 4783, 387, 393,
                                                                       2191, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6007, 0, 6, 4783,
                                                                       1597, 4801, 393, 399,
                                                                       2209, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6043, 0, 6, 4801,
                                                                       1606, 4819, 399, 405,
                                                                       2227, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 6079, 0, 6, 4819,
                                                                       1615, 4837, 405, 411,
                                                                       2245, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 6115, 0, 3, 6,
                                                                       4855, 1633, 4909, 5611,
                                                                       2011, 5647, 423, 441,
                                                                       2263, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 6223, 0, 3, 6,
                                                                       4909, 1660, 4963, 5647,
                                                                       2029, 5683, 441, 459,
                                                                       2317, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 6331, 0, 3, 6,
                                                                       4963, 1687, 5017, 5683,
                                                                       2047, 5719, 459, 477,
                                                                       2371, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 6439, 0, 3, 6,
                                                                       5017, 1714, 5071, 5719,
                                                                       2065, 5755, 477, 495,
                                                                       2425, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 6547, 0, 3, 6,
                                                                       5071, 1741, 5125, 5755,
                                                                       2083, 5791, 495, 513,
                                                                       2479, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 6655, 0, 3, 6,
                                                                       5125, 1768, 5179, 5791,
                                                                       2101, 5827, 513, 531,
                                                                       2533, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 6763, 0, 3, 6,
                                                                       5233, 1822, 5287, 5863,
                                                                       2137, 5899, 567, 585,
                                                                       2587, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 6871, 0, 3, 6,
                                                                       5287, 1849, 5341, 5899,
                                                                       2155, 5935, 585, 603,
                                                                       2641, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 6979, 0, 3, 6,
                                                                       5341, 1876, 5395, 5935,
                                                                       2173, 5971, 603, 621,
                                                                       2695, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 7087, 0, 3, 6,
                                                                       5395, 1903, 5449, 5971,
                                                                       2191, 6007, 621, 639,
                                                                       2749, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 7195, 0, 3, 6,
                                                                       5449, 1930, 5503, 6007,
                                                                       2209, 6043, 639, 657,
                                                                       2803, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 7303, 0, 3, 6,
                                                                       5503, 1957, 5557, 6043,
                                                                       2227, 6079, 657, 675,
                                                                       2857, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7411, 0, 6, 5611,
                                                                       2011, 5647, 711, 721,
                                                                       2911, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7471, 0, 6, 5647,
                                                                       2029, 5683, 721, 731,
                                                                       2941, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7531, 0, 6, 5683,
                                                                       2047, 5719, 731, 741,
                                                                       2971, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7591, 0, 6, 5719,
                                                                       2065, 5755, 741, 751,
                                                                       3001, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7651, 0, 6, 5755,
                                                                       2083, 5791, 751, 761,
                                                                       3031, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7711, 0, 6, 5791,
                                                                       2101, 5827, 761, 771,
                                                                       3061, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7771, 0, 6, 5863,
                                                                       2137, 5899, 791, 801,
                                                                       3091, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7831, 0, 6, 5899,
                                                                       2155, 5935, 801, 811,
                                                                       3121, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7891, 0, 6, 5935,
                                                                       2173, 5971, 811, 821,
                                                                       3151, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 7951, 0, 6, 5971,
                                                                       2191, 6007, 821, 831,
                                                                       3181, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8011, 0, 6, 6007,
                                                                       2209, 6043, 831, 841,
                                                                       3211, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 8071, 0, 6, 6043,
                                                                       2227, 6079, 841, 851,
                                                                       3241, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 8131, 0, 3, 6,
                                                                       4855, 4909, 6115, 2263,
                                                                       6223, 7411, 2911, 7471,
                                                                       871, 901, 3271, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 8311, 0, 3, 6,
                                                                       4909, 4963, 6223, 2317,
                                                                       6331, 7471, 2941, 7531,
                                                                       901, 931, 3361, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 8491, 0, 3, 6,
                                                                       4963, 5017, 6331, 2371,
                                                                       6439, 7531, 2971, 7591,
                                                                       931, 961, 3451, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 8671, 0, 3, 6,
                                                                       5017, 5071, 6439, 2425,
                                                                       6547, 7591, 3001, 7651,
                                                                       961, 991, 3541, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 8851, 0, 3, 6,
                                                                       5071, 5125, 6547, 2479,
                                                                       6655, 7651, 3031, 7711,
                                                                       991, 1021, 3631, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 9031, 0, 3, 6,
                                                                       5233, 5287, 6763, 2587,
                                                                       6871, 7771, 3091, 7831,
                                                                       1081, 1111, 3721, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 9211, 0, 3, 6,
                                                                       5287, 5341, 6871, 2641,
                                                                       6979, 7831, 3121, 7891,
                                                                       1111, 1141, 3811, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 9391, 0, 3, 6,
                                                                       5341, 5395, 6979, 2695,
                                                                       7087, 7891, 3151, 7951,
                                                                       1141, 1171, 3901, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 9571, 0, 3, 6,
                                                                       5395, 5449, 7087, 2749,
                                                                       7195, 7951, 3181, 8011,
                                                                       1171, 1201, 3991, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 9751, 0, 3, 6,
                                                                       5449, 5503, 7195, 2803,
                                                                       7303, 8011, 3211, 8071,
                                                                       1201, 1231, 4081, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9931, 6, 1291,
                                                                       1294, 4183, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9941, 6, 1294,
                                                                       1297, 4189, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9951, 6, 1297,
                                                                       1300, 4195, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9961, 6, 1300,
                                                                       1303, 4201, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9971, 6, 1303,
                                                                       1306, 4207, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9981, 6, 1306,
                                                                       1309, 4213, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9991, 6, 1309,
                                                                       1312, 4219, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10001, 6, 1318,
                                                                       1321, 4237, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10011, 6, 1321,
                                                                       1324, 4243, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10021, 6, 1324,
                                                                       1327, 4249, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10031, 6, 1327,
                                                                       1330, 4255, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10041, 6, 1330,
                                                                       1333, 4261, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10051, 6, 1333,
                                                                       1336, 4267, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 10061, 6, 1336,
                                                                       1339, 4273, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10071, 3, 6, 9931,
                                                                       4183, 9941, 4315, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10101, 3, 6, 9941,
                                                                       4189, 9951, 4333, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10131, 3, 6, 9951,
                                                                       4195, 9961, 4351, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10161, 3, 6, 9961,
                                                                       4201, 9971, 4369, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10191, 3, 6, 9971,
                                                                       4207, 9981, 4387, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10221, 3, 6, 9981,
                                                                       4213, 9991, 4405, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10251, 3, 6,
                                                                       10001, 4237, 10011, 4459,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10281, 3, 6,
                                                                       10011, 4243, 10021, 4477,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10311, 3, 6,
                                                                       10021, 4249, 10031, 4495,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10341, 3, 6,
                                                                       10031, 4255, 10041, 4513,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10371, 3, 6,
                                                                       10041, 4261, 10051, 4531,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 10401, 3, 6,
                                                                       10051, 4267, 10061, 4549,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10431, 0, 6, 9931,
                                                                       4183, 9941, 1489, 1498,
                                                                       4603, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10461, 0, 6, 9941,
                                                                       4189, 9951, 1498, 1507,
                                                                       4621, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10491, 0, 6, 9951,
                                                                       4195, 9961, 1507, 1516,
                                                                       4639, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10521, 0, 6, 9961,
                                                                       4201, 9971, 1516, 1525,
                                                                       4657, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10551, 0, 6, 9971,
                                                                       4207, 9981, 1525, 1534,
                                                                       4675, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10581, 0, 6, 9981,
                                                                       4213, 9991, 1534, 1543,
                                                                       4693, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10611, 0, 6,
                                                                       10001, 4237, 10011, 1561,
                                                                       1570, 4747, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10641, 0, 6,
                                                                       10011, 4243, 10021, 1570,
                                                                       1579, 4765, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10671, 0, 6,
                                                                       10021, 4249, 10031, 1579,
                                                                       1588, 4783, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10701, 0, 6,
                                                                       10031, 4255, 10041, 1588,
                                                                       1597, 4801, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10731, 0, 6,
                                                                       10041, 4261, 10051, 1597,
                                                                       1606, 4819, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 10761, 0, 6,
                                                                       10051, 4267, 10061, 1606,
                                                                       1615, 4837, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10791, 0, 3, 6,
                                                                       10071, 4315, 10101, 10431,
                                                                       4603, 10461, 1633, 1660,
                                                                       4963, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10881, 0, 3, 6,
                                                                       10101, 4333, 10131, 10461,
                                                                       4621, 10491, 1660, 1687,
                                                                       5017, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 10971, 0, 3, 6,
                                                                       10131, 4351, 10161, 10491,
                                                                       4639, 10521, 1687, 1714,
                                                                       5071, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11061, 0, 3, 6,
                                                                       10161, 4369, 10191, 10521,
                                                                       4657, 10551, 1714, 1741,
                                                                       5125, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11151, 0, 3, 6,
                                                                       10191, 4387, 10221, 10551,
                                                                       4675, 10581, 1741, 1768,
                                                                       5179, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11241, 0, 3, 6,
                                                                       10251, 4459, 10281, 10611,
                                                                       4747, 10641, 1822, 1849,
                                                                       5341, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11331, 0, 3, 6,
                                                                       10281, 4477, 10311, 10641,
                                                                       4765, 10671, 1849, 1876,
                                                                       5395, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11421, 0, 3, 6,
                                                                       10311, 4495, 10341, 10671,
                                                                       4783, 10701, 1876, 1903,
                                                                       5449, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11511, 0, 3, 6,
                                                                       10341, 4513, 10371, 10701,
                                                                       4801, 10731, 1903, 1930,
                                                                       5503, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 11601, 0, 3, 6,
                                                                       10371, 4531, 10401, 10731,
                                                                       4819, 10761, 1930, 1957,
                                                                       5557, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 11691, 0, 6,
                                                                       10431, 4603, 10461, 2011,
                                                                       2029, 5683, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 11751, 0, 6,
                                                                       10461, 4621, 10491, 2029,
                                                                       2047, 5719, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 11811, 0, 6,
                                                                       10491, 4639, 10521, 2047,
                                                                       2065, 5755, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 11871, 0, 6,
                                                                       10521, 4657, 10551, 2065,
                                                                       2083, 5791, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 11931, 0, 6,
                                                                       10551, 4675, 10581, 2083,
                                                                       2101, 5827, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 11991, 0, 6,
                                                                       10611, 4747, 10641, 2137,
                                                                       2155, 5935, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12051, 0, 6,
                                                                       10641, 4765, 10671, 2155,
                                                                       2173, 5971, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12111, 0, 6,
                                                                       10671, 4783, 10701, 2173,
                                                                       2191, 6007, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12171, 0, 6,
                                                                       10701, 4801, 10731, 2191,
                                                                       2209, 6043, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 12231, 0, 6,
                                                                       10731, 4819, 10761, 2209,
                                                                       2227, 6079, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 12291, 0, 3, 6,
                                                                       10791, 4963, 10881, 11691,
                                                                       5683, 11751, 2263, 2317,
                                                                       6331, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 12471, 0, 3, 6,
                                                                       10881, 5017, 10971, 11751,
                                                                       5719, 11811, 2317, 2371,
                                                                       6439, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 12651, 0, 3, 6,
                                                                       10971, 5071, 11061, 11811,
                                                                       5755, 11871, 2371, 2425,
                                                                       6547, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 12831, 0, 3, 6,
                                                                       11061, 5125, 11151, 11871,
                                                                       5791, 11931, 2425, 2479,
                                                                       6655, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 13011, 0, 3, 6,
                                                                       11241, 5341, 11331, 11991,
                                                                       5935, 12051, 2587, 2641,
                                                                       6979, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 13191, 0, 3, 6,
                                                                       11331, 5395, 11421, 12051,
                                                                       5971, 12111, 2641, 2695,
                                                                       7087, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 13371, 0, 3, 6,
                                                                       11421, 5449, 11511, 12111,
                                                                       6007, 12171, 2695, 2749,
                                                                       7195, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 13551, 0, 3, 6,
                                                                       11511, 5503, 11601, 12171,
                                                                       6043, 12231, 2749, 2803,
                                                                       7303, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13731, 0, 6,
                                                                       11691, 5683, 11751, 2911,
                                                                       2941, 7531, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13831, 0, 6,
                                                                       11751, 5719, 11811, 2941,
                                                                       2971, 7591, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 13931, 0, 6,
                                                                       11811, 5755, 11871, 2971,
                                                                       3001, 7651, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14031, 0, 6,
                                                                       11871, 5791, 11931, 3001,
                                                                       3031, 7711, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14131, 0, 6,
                                                                       11991, 5935, 12051, 3091,
                                                                       3121, 7891, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14231, 0, 6,
                                                                       12051, 5971, 12111, 3121,
                                                                       3151, 7951, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14331, 0, 6,
                                                                       12111, 6007, 12171, 3151,
                                                                       3181, 8011, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 14431, 0, 6,
                                                                       12171, 6043, 12231, 3181,
                                                                       3211, 8071, ncols, gamma,
                                                                       p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 14531, 0, 3, 6,
                                                                       10791, 10881, 12291, 6331,
                                                                       12471, 13731, 7531, 13831,
                                                                       3271, 3361, 8491, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 14831, 0, 3, 6,
                                                                       10881, 10971, 12471, 6439,
                                                                       12651, 13831, 7591, 13931,
                                                                       3361, 3451, 8671, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 15131, 0, 3, 6,
                                                                       10971, 11061, 12651, 6547,
                                                                       12831, 13931, 7651, 14031,
                                                                       3451, 3541, 8851, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 15431, 0, 3, 6,
                                                                       11241, 11331, 13011, 6979,
                                                                       13191, 14131, 7891, 14231,
                                                                       3721, 3811, 9391, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 15731, 0, 3, 6,
                                                                       11331, 11421, 13191, 7087,
                                                                       13371, 14231, 7951, 14331,
                                                                       3811, 3901, 9571, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 16031, 0, 3, 6,
                                                                       11421, 11511, 13371, 7195,
                                                                       13551, 14331, 8011, 14431,
                                                                       3901, 3991, 9751, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16331, 6, 4171,
                                                                       4177, 9931, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16346, 6, 4177,
                                                                       4183, 9941, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16361, 6, 4183,
                                                                       4189, 9951, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16376, 6, 4189,
                                                                       4195, 9961, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16391, 6, 4195,
                                                                       4201, 9971, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16406, 6, 4201,
                                                                       4207, 9981, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16421, 6, 4207,
                                                                       4213, 9991, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16436, 6, 4225,
                                                                       4231, 10001, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16451, 6, 4231,
                                                                       4237, 10011, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16466, 6, 4237,
                                                                       4243, 10021, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16481, 6, 4243,
                                                                       4249, 10031, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16496, 6, 4249,
                                                                       4255, 10041, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16511, 6, 4255,
                                                                       4261, 10051, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 16526, 6, 4261,
                                                                       4267, 10061, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16541, 3, 6,
                                                                       16331, 9931, 16346, 4279,
                                                                       4297, 10071, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16586, 3, 6,
                                                                       16346, 9941, 16361, 4297,
                                                                       4315, 10101, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16631, 3, 6,
                                                                       16361, 9951, 16376, 4315,
                                                                       4333, 10131, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16676, 3, 6,
                                                                       16376, 9961, 16391, 4333,
                                                                       4351, 10161, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16721, 3, 6,
                                                                       16391, 9971, 16406, 4351,
                                                                       4369, 10191, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16766, 3, 6,
                                                                       16406, 9981, 16421, 4369,
                                                                       4387, 10221, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16811, 3, 6,
                                                                       16436, 10001, 16451, 4423,
                                                                       4441, 10251, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16856, 3, 6,
                                                                       16451, 10011, 16466, 4441,
                                                                       4459, 10281, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16901, 3, 6,
                                                                       16466, 10021, 16481, 4459,
                                                                       4477, 10311, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16946, 3, 6,
                                                                       16481, 10031, 16496, 4477,
                                                                       4495, 10341, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 16991, 3, 6,
                                                                       16496, 10041, 16511, 4495,
                                                                       4513, 10371, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 17036, 3, 6,
                                                                       16511, 10051, 16526, 4513,
                                                                       4531, 10401, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17081, 0, 6,
                                                                       16331, 9931, 16346, 4567,
                                                                       4585, 10431, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17126, 0, 6,
                                                                       16346, 9941, 16361, 4585,
                                                                       4603, 10461, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17171, 0, 6,
                                                                       16361, 9951, 16376, 4603,
                                                                       4621, 10491, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17216, 0, 6,
                                                                       16376, 9961, 16391, 4621,
                                                                       4639, 10521, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17261, 0, 6,
                                                                       16391, 9971, 16406, 4639,
                                                                       4657, 10551, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17306, 0, 6,
                                                                       16406, 9981, 16421, 4657,
                                                                       4675, 10581, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17351, 0, 6,
                                                                       16436, 10001, 16451, 4711,
                                                                       4729, 10611, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17396, 0, 6,
                                                                       16451, 10011, 16466, 4729,
                                                                       4747, 10641, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17441, 0, 6,
                                                                       16466, 10021, 16481, 4747,
                                                                       4765, 10671, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17486, 0, 6,
                                                                       16481, 10031, 16496, 4765,
                                                                       4783, 10701, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17531, 0, 6,
                                                                       16496, 10041, 16511, 4783,
                                                                       4801, 10731, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 17576, 0, 6,
                                                                       16511, 10051, 16526, 4801,
                                                                       4819, 10761, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 17621, 0, 3, 6,
                                                                       16541, 10071, 16586,
                                                                       17081, 10431, 17126, 4855,
                                                                       4909, 10791, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 17756, 0, 3, 6,
                                                                       16586, 10101, 16631,
                                                                       17126, 10461, 17171, 4909,
                                                                       4963, 10881, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 17891, 0, 3, 6,
                                                                       16631, 10131, 16676,
                                                                       17171, 10491, 17216, 4963,
                                                                       5017, 10971, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 18026, 0, 3, 6,
                                                                       16676, 10161, 16721,
                                                                       17216, 10521, 17261, 5017,
                                                                       5071, 11061, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 18161, 0, 3, 6,
                                                                       16721, 10191, 16766,
                                                                       17261, 10551, 17306, 5071,
                                                                       5125, 11151, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 18296, 0, 3, 6,
                                                                       16811, 10251, 16856,
                                                                       17351, 10611, 17396, 5233,
                                                                       5287, 11241, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 18431, 0, 3, 6,
                                                                       16856, 10281, 16901,
                                                                       17396, 10641, 17441, 5287,
                                                                       5341, 11331, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 18566, 0, 3, 6,
                                                                       16901, 10311, 16946,
                                                                       17441, 10671, 17486, 5341,
                                                                       5395, 11421, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 18701, 0, 3, 6,
                                                                       16946, 10341, 16991,
                                                                       17486, 10701, 17531, 5395,
                                                                       5449, 11511, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 18836, 0, 3, 6,
                                                                       16991, 10371, 17036,
                                                                       17531, 10731, 17576, 5449,
                                                                       5503, 11601, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 18971, 0, 6,
                                                                       17081, 10431, 17126, 5611,
                                                                       5647, 11691, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 19061, 0, 6,
                                                                       17126, 10461, 17171, 5647,
                                                                       5683, 11751, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 19151, 0, 6,
                                                                       17171, 10491, 17216, 5683,
                                                                       5719, 11811, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 19241, 0, 6,
                                                                       17216, 10521, 17261, 5719,
                                                                       5755, 11871, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 19331, 0, 6,
                                                                       17261, 10551, 17306, 5755,
                                                                       5791, 11931, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 19421, 0, 6,
                                                                       17351, 10611, 17396, 5863,
                                                                       5899, 11991, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 19511, 0, 6,
                                                                       17396, 10641, 17441, 5899,
                                                                       5935, 12051, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 19601, 0, 6,
                                                                       17441, 10671, 17486, 5935,
                                                                       5971, 12111, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 19691, 0, 6,
                                                                       17486, 10701, 17531, 5971,
                                                                       6007, 12171, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 19781, 0, 6,
                                                                       17531, 10731, 17576, 6007,
                                                                       6043, 12231, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 19871, 0, 3, 6,
                                                                       17621, 10791, 17756,
                                                                       18971, 11691, 19061, 6115,
                                                                       6223, 12291, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 20141, 0, 3, 6,
                                                                       17756, 10881, 17891,
                                                                       19061, 11751, 19151, 6223,
                                                                       6331, 12471, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 20411, 0, 3, 6,
                                                                       17891, 10971, 18026,
                                                                       19151, 11811, 19241, 6331,
                                                                       6439, 12651, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 20681, 0, 3, 6,
                                                                       18026, 11061, 18161,
                                                                       19241, 11871, 19331, 6439,
                                                                       6547, 12831, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 20951, 0, 3, 6,
                                                                       18296, 11241, 18431,
                                                                       19421, 11991, 19511, 6763,
                                                                       6871, 13011, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 21221, 0, 3, 6,
                                                                       18431, 11331, 18566,
                                                                       19511, 12051, 19601, 6871,
                                                                       6979, 13191, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 21491, 0, 3, 6,
                                                                       18566, 11421, 18701,
                                                                       19601, 12111, 19691, 6979,
                                                                       7087, 13371, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 21761, 0, 3, 6,
                                                                       18701, 11511, 18836,
                                                                       19691, 12171, 19781, 7087,
                                                                       7195, 13551, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22031, 0, 6,
                                                                       18971, 11691, 19061, 7411,
                                                                       7471, 13731, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22181, 0, 6,
                                                                       19061, 11751, 19151, 7471,
                                                                       7531, 13831, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22331, 0, 6,
                                                                       19151, 11811, 19241, 7531,
                                                                       7591, 13931, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22481, 0, 6,
                                                                       19241, 11871, 19331, 7591,
                                                                       7651, 14031, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22631, 0, 6,
                                                                       19421, 11991, 19511, 7771,
                                                                       7831, 14131, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22781, 0, 6,
                                                                       19511, 12051, 19601, 7831,
                                                                       7891, 14231, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 22931, 0, 6,
                                                                       19601, 12111, 19691, 7891,
                                                                       7951, 14331, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 23081, 0, 6,
                                                                       19691, 12171, 19781, 7951,
                                                                       8011, 14431, ncols, gamma,
                                                                       p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 23231, 0, 3, 6,
                                                                       17621, 17756, 19871,
                                                                       12291, 20141, 22031,
                                                                       13731, 22181, 8131, 8311,
                                                                       14531, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 23681, 0, 3, 6,
                                                                       17756, 17891, 20141,
                                                                       12471, 20411, 22181,
                                                                       13831, 22331, 8311, 8491,
                                                                       14831, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 24131, 0, 3, 6,
                                                                       17891, 18026, 20411,
                                                                       12651, 20681, 22331,
                                                                       13931, 22481, 8491, 8671,
                                                                       15131, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 24581, 0, 3, 6,
                                                                       18296, 18431, 20951,
                                                                       13011, 21221, 22631,
                                                                       14131, 22781, 9031, 9211,
                                                                       15431, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 25031, 0, 3, 6,
                                                                       18431, 18566, 21221,
                                                                       13191, 21491, 22781,
                                                                       14231, 22931, 9211, 9391,
                                                                       15731, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 25481, 0, 3, 6,
                                                                       18566, 18701, 21491,
                                                                       13371, 21761, 22931,
                                                                       14331, 23081, 9391, 9571,
                                                                       16031, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25931, 6, 9931,
                                                                       9941, 16361, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25952, 6, 9941,
                                                                       9951, 16376, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25973, 6, 9951,
                                                                       9961, 16391, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 25994, 6, 9961,
                                                                       9971, 16406, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26015, 6, 9971,
                                                                       9981, 16421, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26036, 6, 10001,
                                                                       10011, 16466, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26057, 6, 10011,
                                                                       10021, 16481, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26078, 6, 10021,
                                                                       10031, 16496, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26099, 6, 10031,
                                                                       10041, 16511, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 26120, 6, 10041,
                                                                       10051, 16526, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26141, 3, 6,
                                                                       25931, 16361, 25952,
                                                                       10071, 10101, 16631,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26204, 3, 6,
                                                                       25952, 16376, 25973,
                                                                       10101, 10131, 16676,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26267, 3, 6,
                                                                       25973, 16391, 25994,
                                                                       10131, 10161, 16721,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26330, 3, 6,
                                                                       25994, 16406, 26015,
                                                                       10161, 10191, 16766,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26393, 3, 6,
                                                                       26036, 16466, 26057,
                                                                       10251, 10281, 16901,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26456, 3, 6,
                                                                       26057, 16481, 26078,
                                                                       10281, 10311, 16946,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26519, 3, 6,
                                                                       26078, 16496, 26099,
                                                                       10311, 10341, 16991,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 26582, 3, 6,
                                                                       26099, 16511, 26120,
                                                                       10341, 10371, 17036,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26645, 0, 6,
                                                                       25931, 16361, 25952,
                                                                       10431, 10461, 17171,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26708, 0, 6,
                                                                       25952, 16376, 25973,
                                                                       10461, 10491, 17216,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26771, 0, 6,
                                                                       25973, 16391, 25994,
                                                                       10491, 10521, 17261,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26834, 0, 6,
                                                                       25994, 16406, 26015,
                                                                       10521, 10551, 17306,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26897, 0, 6,
                                                                       26036, 16466, 26057,
                                                                       10611, 10641, 17441,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 26960, 0, 6,
                                                                       26057, 16481, 26078,
                                                                       10641, 10671, 17486,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 27023, 0, 6,
                                                                       26078, 16496, 26099,
                                                                       10671, 10701, 17531,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 27086, 0, 6,
                                                                       26099, 16511, 26120,
                                                                       10701, 10731, 17576,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 27149, 0, 3, 6,
                                                                       26141, 16631, 26204,
                                                                       26645, 17171, 26708,
                                                                       10791, 10881, 17891,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 27338, 0, 3, 6,
                                                                       26204, 16676, 26267,
                                                                       26708, 17216, 26771,
                                                                       10881, 10971, 18026,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 27527, 0, 3, 6,
                                                                       26267, 16721, 26330,
                                                                       26771, 17261, 26834,
                                                                       10971, 11061, 18161,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 27716, 0, 3, 6,
                                                                       26393, 16901, 26456,
                                                                       26897, 17441, 26960,
                                                                       11241, 11331, 18566,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 27905, 0, 3, 6,
                                                                       26456, 16946, 26519,
                                                                       26960, 17486, 27023,
                                                                       11331, 11421, 18701,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 28094, 0, 3, 6,
                                                                       26519, 16991, 26582,
                                                                       27023, 17531, 27086,
                                                                       11421, 11511, 18836,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 28283, 0, 6,
                                                                       26645, 17171, 26708,
                                                                       11691, 11751, 19151,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 28409, 0, 6,
                                                                       26708, 17216, 26771,
                                                                       11751, 11811, 19241,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 28535, 0, 6,
                                                                       26771, 17261, 26834,
                                                                       11811, 11871, 19331,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 28661, 0, 6,
                                                                       26897, 17441, 26960,
                                                                       11991, 12051, 19601,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 28787, 0, 6,
                                                                       26960, 17486, 27023,
                                                                       12051, 12111, 19691,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 28913, 0, 6,
                                                                       27023, 17531, 27086,
                                                                       12111, 12171, 19781,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 29039, 0, 3, 6,
                                                                       27149, 17891, 27338,
                                                                       28283, 19151, 28409,
                                                                       12291, 12471, 20411,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 29417, 0, 3, 6,
                                                                       27338, 18026, 27527,
                                                                       28409, 19241, 28535,
                                                                       12471, 12651, 20681,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 29795, 0, 3, 6,
                                                                       27716, 18566, 27905,
                                                                       28661, 19601, 28787,
                                                                       13011, 13191, 21491,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 30173, 0, 3, 6,
                                                                       27905, 18701, 28094,
                                                                       28787, 19691, 28913,
                                                                       13191, 13371, 21761,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 30551, 0, 6,
                                                                       28283, 19151, 28409,
                                                                       13731, 13831, 22331,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 30761, 0, 6,
                                                                       28409, 19241, 28535,
                                                                       13831, 13931, 22481,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 30971, 0, 6,
                                                                       28661, 19601, 28787,
                                                                       14131, 14231, 22931,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 31181, 0, 6,
                                                                       28787, 19691, 28913,
                                                                       14231, 14331, 23081,
                                                                       ncols, gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 31391, 0, 3, 6,
                                                                       27149, 27338, 29039,
                                                                       20411, 29417, 30551,
                                                                       22331, 30761, 14531,
                                                                       14831, 24131, ncols,
                                                                       gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 32021, 0, 3, 6,
                                                                       27716, 27905, 29795,
                                                                       21491, 30173, 30971,
                                                                       22931, 31181, 15431,
                                                                       15731, 25481, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32651, 6, 16331,
                                                                       16346, 25931, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32679, 6, 16346,
                                                                       16361, 25952, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32707, 6, 16361,
                                                                       16376, 25973, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32735, 6, 16376,
                                                                       16391, 25994, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32763, 6, 16391,
                                                                       16406, 26015, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32791, 6, 16436,
                                                                       16451, 26036, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32819, 6, 16451,
                                                                       16466, 26057, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32847, 6, 16466,
                                                                       16481, 26078, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32875, 6, 16481,
                                                                       16496, 26099, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 32903, 6, 16496,
                                                                       16511, 26120, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 32931, 3, 6,
                                                                       32651, 25931, 32679,
                                                                       16541, 16586, 26141,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33015, 3, 6,
                                                                       32679, 25952, 32707,
                                                                       16586, 16631, 26204,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33099, 3, 6,
                                                                       32707, 25973, 32735,
                                                                       16631, 16676, 26267,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33183, 3, 6,
                                                                       32735, 25994, 32763,
                                                                       16676, 16721, 26330,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33267, 3, 6,
                                                                       32791, 26036, 32819,
                                                                       16811, 16856, 26393,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33351, 3, 6,
                                                                       32819, 26057, 32847,
                                                                       16856, 16901, 26456,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33435, 3, 6,
                                                                       32847, 26078, 32875,
                                                                       16901, 16946, 26519,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 33519, 3, 6,
                                                                       32875, 26099, 32903,
                                                                       16946, 16991, 26582,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 33603, 0, 6,
                                                                       32651, 25931, 32679,
                                                                       17081, 17126, 26645,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 33687, 0, 6,
                                                                       32679, 25952, 32707,
                                                                       17126, 17171, 26708,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 33771, 0, 6,
                                                                       32707, 25973, 32735,
                                                                       17171, 17216, 26771,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 33855, 0, 6,
                                                                       32735, 25994, 32763,
                                                                       17216, 17261, 26834,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 33939, 0, 6,
                                                                       32791, 26036, 32819,
                                                                       17351, 17396, 26897,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 34023, 0, 6,
                                                                       32819, 26057, 32847,
                                                                       17396, 17441, 26960,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 34107, 0, 6,
                                                                       32847, 26078, 32875,
                                                                       17441, 17486, 27023,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 34191, 0, 6,
                                                                       32875, 26099, 32903,
                                                                       17486, 17531, 27086,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 34275, 0, 3, 6,
                                                                       32931, 26141, 33015,
                                                                       33603, 26645, 33687,
                                                                       17621, 17756, 27149,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 34527, 0, 3, 6,
                                                                       33015, 26204, 33099,
                                                                       33687, 26708, 33771,
                                                                       17756, 17891, 27338,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 34779, 0, 3, 6,
                                                                       33099, 26267, 33183,
                                                                       33771, 26771, 33855,
                                                                       17891, 18026, 27527,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 35031, 0, 3, 6,
                                                                       33267, 26393, 33351,
                                                                       33939, 26897, 34023,
                                                                       18296, 18431, 27716,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 35283, 0, 3, 6,
                                                                       33351, 26456, 33435,
                                                                       34023, 26960, 34107,
                                                                       18431, 18566, 27905,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 35535, 0, 3, 6,
                                                                       33435, 26519, 33519,
                                                                       34107, 27023, 34191,
                                                                       18566, 18701, 28094,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 35787, 0, 6,
                                                                       33603, 26645, 33687,
                                                                       18971, 19061, 28283,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 35955, 0, 6,
                                                                       33687, 26708, 33771,
                                                                       19061, 19151, 28409,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 36123, 0, 6,
                                                                       33771, 26771, 33855,
                                                                       19151, 19241, 28535,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 36291, 0, 6,
                                                                       33939, 26897, 34023,
                                                                       19421, 19511, 28661,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 36459, 0, 6,
                                                                       34023, 26960, 34107,
                                                                       19511, 19601, 28787,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 36627, 0, 6,
                                                                       34107, 27023, 34191,
                                                                       19601, 19691, 28913,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 36795, 0, 3, 6,
                                                                       34275, 27149, 34527,
                                                                       35787, 28283, 35955,
                                                                       19871, 20141, 29039,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 37299, 0, 3, 6,
                                                                       34527, 27338, 34779,
                                                                       35955, 28409, 36123,
                                                                       20141, 20411, 29417,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 37803, 0, 3, 6,
                                                                       35031, 27716, 35283,
                                                                       36291, 28661, 36459,
                                                                       20951, 21221, 29795,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 38307, 0, 3, 6,
                                                                       35283, 27905, 35535,
                                                                       36459, 28787, 36627,
                                                                       21221, 21491, 30173,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 38811, 0, 6,
                                                                       35787, 28283, 35955,
                                                                       22031, 22181, 30551,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 39091, 0, 6,
                                                                       35955, 28409, 36123,
                                                                       22181, 22331, 30761,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 39371, 0, 6,
                                                                       36291, 28661, 36459,
                                                                       22631, 22781, 30971,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 39651, 0, 6,
                                                                       36459, 28787, 36627,
                                                                       22781, 22931, 31181,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpi_three_center_electron_repulsion_0(buffer, 39931, 0, 3, 6,
                                                                       34275, 34527, 36795,
                                                                       29039, 37299, 38811,
                                                                       30551, 39091, 23231,
                                                                       23681, 31391, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpi_three_center_electron_repulsion_0(buffer, 40771, 0, 3, 6,
                                                                       35031, 35283, 37803,
                                                                       29795, 38307, 39371,
                                                                       30971, 39651, 24581,
                                                                       25031, 32021, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 41611, 40771, 10, 28, ncols, beta);

                    simdgeo::geom_s_y(buffer, 41891, 40771, 10, 28, ncols, beta);

                    simdgeo::geom_s_z(buffer, 42171, 40771, 10, 28, ncols, beta);

                    simdgeo::geom_s_x(buffer, 42451, 39931, 10, 28, ncols, beta);

                    simdgeo::geom_s_y(buffer, 42731, 39931, 10, 28, ncols, beta);

                    simdgeo::geom_s_z(buffer, 43011, 39931, 10, 28, ncols, beta);

                    simdfunc::contract_primitives(buffer, 43291, 41611, 1680, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 44971, 43291, 10, 1, nmax);

        simdtrf::transform_f_outer(values + n * npairs, nvalues, buffer, 44971, 13, nmax);

        simdtrf::transform_i_inner(buffer, 44971, 43571, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 91 * nvalues + n * npairs, nvalues, buffer, 44971,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 44971, 43851, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 182 * nvalues + n * npairs, nvalues, buffer, 44971,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 44971, 44131, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 273 * nvalues + n * npairs, nvalues, buffer, 44971,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 44971, 44411, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 364 * nvalues + n * npairs, nvalues, buffer, 44971,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 44971, 44691, 10, 1, nmax);

        simdtrf::transform_f_outer(values + 455 * nvalues + n * npairs, nvalues, buffer, 44971,
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
