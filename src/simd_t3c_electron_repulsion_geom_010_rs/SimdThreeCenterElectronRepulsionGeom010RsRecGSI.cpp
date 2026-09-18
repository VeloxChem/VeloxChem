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


#include "SimdThreeCenterElectronRepulsionGeom010RsRecGSI.hpp"

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
#include "SimdThreeCenterElectronRepulsionVrrRecGPD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGPS.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSD.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSF.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSG.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSH.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSI.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSP.hpp"
#include "SimdThreeCenterElectronRepulsionVrrRecGSS.hpp"
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
#include "SimdTransformG.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_010_gsi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_010_gsi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 83926, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 702 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 83926, 81211, 2520, dimensions);

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

                    simdfunc::compute_full_t3c_erf_boys_function(buffer, coordinates, 9, 6, 11,
                                                                 ncols, fj, i * nprim_b + j, fq,
                                                                 omega);

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 22, 6, 11,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 35, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 38, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 41, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 44, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 47, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 50, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 53, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 56, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 59, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 62, 3, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 65, 3, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 68, 3, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 71, 3, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 74, 3, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 77, 3, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 80, 3, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 83, 3, 6, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 86, 3, 6, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 89, 3, 6, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 92, 3, 6, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 95, 3, 6, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 98, 3, 6, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 101, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 104, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 107, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 110, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 113, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 116, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 119, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 122, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 125, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 128, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 131, 0, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 134, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 137, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 140, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 143, 0, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 146, 0, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 149, 0, 6, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 152, 0, 6, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 155, 0, 6, 30, 31,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 158, 0, 6, 31, 32,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 161, 0, 6, 32, 33,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 164, 0, 6, 33, 34,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 167, 0, 6, 10, 11,
                                                                       35, 38, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 176, 0, 6, 11, 12,
                                                                       38, 41, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 185, 0, 6, 12, 13,
                                                                       41, 44, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 194, 0, 6, 13, 14,
                                                                       44, 47, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 203, 0, 6, 14, 15,
                                                                       47, 50, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 212, 0, 6, 15, 16,
                                                                       50, 53, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 221, 0, 6, 16, 17,
                                                                       53, 56, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 230, 0, 6, 17, 18,
                                                                       56, 59, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 239, 0, 6, 18, 19,
                                                                       59, 62, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 248, 0, 6, 19, 20,
                                                                       62, 65, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 257, 0, 6, 23, 24,
                                                                       68, 71, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 266, 0, 6, 24, 25,
                                                                       71, 74, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 275, 0, 6, 25, 26,
                                                                       74, 77, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 284, 0, 6, 26, 27,
                                                                       77, 80, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 293, 0, 6, 27, 28,
                                                                       80, 83, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 302, 0, 6, 28, 29,
                                                                       83, 86, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 311, 0, 6, 29, 30,
                                                                       86, 89, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 320, 0, 6, 30, 31,
                                                                       89, 92, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 329, 0, 6, 31, 32,
                                                                       92, 95, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 338, 0, 6, 32, 33,
                                                                       95, 98, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 347, 0, 6, 10, 11,
                                                                       101, 104, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 353, 0, 6, 11, 12,
                                                                       104, 107, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 359, 0, 6, 12, 13,
                                                                       107, 110, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 365, 0, 6, 13, 14,
                                                                       110, 113, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 371, 0, 6, 14, 15,
                                                                       113, 116, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 377, 0, 6, 15, 16,
                                                                       116, 119, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 383, 0, 6, 16, 17,
                                                                       119, 122, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 389, 0, 6, 17, 18,
                                                                       122, 125, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 395, 0, 6, 18, 19,
                                                                       125, 128, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 401, 0, 6, 19, 20,
                                                                       128, 131, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 407, 0, 6, 23, 24,
                                                                       134, 137, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 413, 0, 6, 24, 25,
                                                                       137, 140, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 419, 0, 6, 25, 26,
                                                                       140, 143, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 425, 0, 6, 26, 27,
                                                                       143, 146, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 431, 0, 6, 27, 28,
                                                                       146, 149, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 437, 0, 6, 28, 29,
                                                                       149, 152, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 443, 0, 6, 29, 30,
                                                                       152, 155, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 449, 0, 6, 30, 31,
                                                                       155, 158, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 455, 0, 6, 31, 32,
                                                                       158, 161, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 461, 0, 6, 32, 33,
                                                                       161, 164, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 467, 0, 3, 6, 101,
                                                                       104, 167, 176, 347, 353,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 485, 0, 3, 6, 104,
                                                                       107, 176, 185, 353, 359,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 503, 0, 3, 6, 107,
                                                                       110, 185, 194, 359, 365,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 521, 0, 3, 6, 110,
                                                                       113, 194, 203, 365, 371,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 539, 0, 3, 6, 113,
                                                                       116, 203, 212, 371, 377,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 557, 0, 3, 6, 116,
                                                                       119, 212, 221, 377, 383,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 575, 0, 3, 6, 119,
                                                                       122, 221, 230, 383, 389,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 593, 0, 3, 6, 122,
                                                                       125, 230, 239, 389, 395,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 611, 0, 3, 6, 125,
                                                                       128, 239, 248, 395, 401,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 629, 0, 3, 6, 134,
                                                                       137, 257, 266, 407, 413,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 647, 0, 3, 6, 137,
                                                                       140, 266, 275, 413, 419,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 665, 0, 3, 6, 140,
                                                                       143, 275, 284, 419, 425,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 683, 0, 3, 6, 143,
                                                                       146, 284, 293, 425, 431,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 701, 0, 3, 6, 146,
                                                                       149, 293, 302, 431, 437,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 719, 0, 3, 6, 149,
                                                                       152, 302, 311, 437, 443,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 737, 0, 3, 6, 152,
                                                                       155, 311, 320, 443, 449,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 755, 0, 3, 6, 155,
                                                                       158, 320, 329, 449, 455,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 773, 0, 3, 6, 158,
                                                                       161, 329, 338, 455, 461,
                                                                       ncols, gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 791, 0, 6, 101,
                                                                       104, 347, 353, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 801, 0, 6, 104,
                                                                       107, 353, 359, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 811, 0, 6, 107,
                                                                       110, 359, 365, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 821, 0, 6, 110,
                                                                       113, 365, 371, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 831, 0, 6, 113,
                                                                       116, 371, 377, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 841, 0, 6, 116,
                                                                       119, 377, 383, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 851, 0, 6, 119,
                                                                       122, 383, 389, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 861, 0, 6, 122,
                                                                       125, 389, 395, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 871, 0, 6, 125,
                                                                       128, 395, 401, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 881, 0, 6, 134,
                                                                       137, 407, 413, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 891, 0, 6, 137,
                                                                       140, 413, 419, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 901, 0, 6, 140,
                                                                       143, 419, 425, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 911, 0, 6, 143,
                                                                       146, 425, 431, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 921, 0, 6, 146,
                                                                       149, 431, 437, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 931, 0, 6, 149,
                                                                       152, 437, 443, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 941, 0, 6, 152,
                                                                       155, 443, 449, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 951, 0, 6, 155,
                                                                       158, 449, 455, ncols,
                                                                       gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 961, 0, 6, 158,
                                                                       161, 455, 461, ncols,
                                                                       gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 971, 0, 3, 6, 347,
                                                                       353, 467, 485, 791, 801,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1001, 0, 3, 6,
                                                                       353, 359, 485, 503, 801,
                                                                       811, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1031, 0, 3, 6,
                                                                       359, 365, 503, 521, 811,
                                                                       821, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1061, 0, 3, 6,
                                                                       365, 371, 521, 539, 821,
                                                                       831, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1091, 0, 3, 6,
                                                                       371, 377, 539, 557, 831,
                                                                       841, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1121, 0, 3, 6,
                                                                       377, 383, 557, 575, 841,
                                                                       851, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1151, 0, 3, 6,
                                                                       383, 389, 575, 593, 851,
                                                                       861, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1181, 0, 3, 6,
                                                                       389, 395, 593, 611, 861,
                                                                       871, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1211, 0, 3, 6,
                                                                       407, 413, 629, 647, 881,
                                                                       891, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1241, 0, 3, 6,
                                                                       413, 419, 647, 665, 891,
                                                                       901, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1271, 0, 3, 6,
                                                                       419, 425, 665, 683, 901,
                                                                       911, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1301, 0, 3, 6,
                                                                       425, 431, 683, 701, 911,
                                                                       921, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1331, 0, 3, 6,
                                                                       431, 437, 701, 719, 921,
                                                                       931, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1361, 0, 3, 6,
                                                                       437, 443, 719, 737, 931,
                                                                       941, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1391, 0, 3, 6,
                                                                       443, 449, 737, 755, 941,
                                                                       951, ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 1421, 0, 3, 6,
                                                                       449, 455, 755, 773, 951,
                                                                       961, ncols, gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1451, 0, 6, 347,
                                                                       353, 791, 801, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1466, 0, 6, 353,
                                                                       359, 801, 811, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1481, 0, 6, 359,
                                                                       365, 811, 821, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1496, 0, 6, 365,
                                                                       371, 821, 831, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1511, 0, 6, 371,
                                                                       377, 831, 841, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1526, 0, 6, 377,
                                                                       383, 841, 851, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1541, 0, 6, 383,
                                                                       389, 851, 861, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1556, 0, 6, 389,
                                                                       395, 861, 871, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1571, 0, 6, 407,
                                                                       413, 881, 891, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1586, 0, 6, 413,
                                                                       419, 891, 901, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1601, 0, 6, 419,
                                                                       425, 901, 911, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1616, 0, 6, 425,
                                                                       431, 911, 921, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1631, 0, 6, 431,
                                                                       437, 921, 931, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1646, 0, 6, 437,
                                                                       443, 931, 941, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1661, 0, 6, 443,
                                                                       449, 941, 951, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 1676, 0, 6, 449,
                                                                       455, 951, 961, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1691, 0, 3, 6,
                                                                       467, 485, 791, 801, 971,
                                                                       1001, 1451, 1466, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1736, 0, 3, 6,
                                                                       485, 503, 801, 811, 1001,
                                                                       1031, 1466, 1481, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1781, 0, 3, 6,
                                                                       503, 521, 811, 821, 1031,
                                                                       1061, 1481, 1496, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1826, 0, 3, 6,
                                                                       521, 539, 821, 831, 1061,
                                                                       1091, 1496, 1511, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1871, 0, 3, 6,
                                                                       539, 557, 831, 841, 1091,
                                                                       1121, 1511, 1526, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1916, 0, 3, 6,
                                                                       557, 575, 841, 851, 1121,
                                                                       1151, 1526, 1541, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1961, 0, 3, 6,
                                                                       575, 593, 851, 861, 1151,
                                                                       1181, 1541, 1556, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 2006, 0, 3, 6,
                                                                       629, 647, 881, 891, 1211,
                                                                       1241, 1571, 1586, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 2051, 0, 3, 6,
                                                                       647, 665, 891, 901, 1241,
                                                                       1271, 1586, 1601, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 2096, 0, 3, 6,
                                                                       665, 683, 901, 911, 1271,
                                                                       1301, 1601, 1616, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 2141, 0, 3, 6,
                                                                       683, 701, 911, 921, 1301,
                                                                       1331, 1616, 1631, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 2186, 0, 3, 6,
                                                                       701, 719, 921, 931, 1331,
                                                                       1361, 1631, 1646, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 2231, 0, 3, 6,
                                                                       719, 737, 931, 941, 1361,
                                                                       1391, 1646, 1661, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 2276, 0, 3, 6,
                                                                       737, 755, 941, 951, 1391,
                                                                       1421, 1661, 1676, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2321, 6, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2324, 6, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2327, 6, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2330, 6, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2333, 6, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2336, 6, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2339, 6, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2342, 6, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2345, 6, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2348, 6, 21,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2351, 6, 25,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2354, 6, 26,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2357, 6, 27,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2360, 6, 28,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2363, 6, 29,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2366, 6, 30,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2369, 6, 31,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2372, 6, 32,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2375, 6, 33,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 2378, 6, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2381, 6, 12, 41,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2390, 6, 13, 44,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2399, 6, 14, 47,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2408, 6, 15, 50,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2417, 6, 16, 53,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2426, 6, 17, 56,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2435, 6, 18, 59,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2444, 6, 19, 62,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2453, 6, 20, 65,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2462, 6, 25, 74,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2471, 6, 26, 77,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2480, 6, 27, 80,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2489, 6, 28, 83,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2498, 6, 29, 86,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2507, 6, 30, 89,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2516, 6, 31, 92,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2525, 6, 32, 95,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 2534, 6, 33, 98,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2543, 6, 12, 107,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2552, 6, 13, 110,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2561, 6, 14, 113,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2570, 6, 15, 116,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2579, 6, 16, 119,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2588, 6, 17, 122,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2597, 6, 18, 125,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2606, 6, 19, 128,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2615, 6, 20, 131,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2624, 6, 25, 140,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2633, 6, 26, 143,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2642, 6, 27, 146,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2651, 6, 28, 149,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2660, 6, 29, 152,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2669, 6, 30, 155,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2678, 6, 31, 158,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2687, 6, 32, 161,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 2696, 6, 33, 164,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2705, 6, 41, 107,
                                                                       185, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2732, 6, 44, 110,
                                                                       194, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2759, 6, 47, 113,
                                                                       203, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2786, 6, 50, 116,
                                                                       212, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2813, 6, 53, 119,
                                                                       221, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2840, 6, 56, 122,
                                                                       230, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2867, 6, 59, 125,
                                                                       239, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2894, 6, 62, 128,
                                                                       248, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2921, 6, 74, 140,
                                                                       275, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2948, 6, 77, 143,
                                                                       284, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 2975, 6, 80, 146,
                                                                       293, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3002, 6, 83, 149,
                                                                       302, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3029, 6, 86, 152,
                                                                       311, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3056, 6, 89, 155,
                                                                       320, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3083, 6, 92, 158,
                                                                       329, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 3110, 6, 95, 161,
                                                                       338, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3137, 6, 107, 359,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3155, 6, 110, 365,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3173, 6, 113, 371,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3191, 6, 116, 377,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3209, 6, 119, 383,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3227, 6, 122, 389,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3245, 6, 125, 395,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3263, 6, 128, 401,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3281, 6, 140, 419,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3299, 6, 143, 425,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3317, 6, 146, 431,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3335, 6, 149, 437,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3353, 6, 152, 443,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3371, 6, 155, 449,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3389, 6, 158, 455,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 3407, 6, 161, 461,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3425, 0, 6, 2705,
                                                                       185, 2732, 359, 503,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3479, 0, 6, 2732,
                                                                       194, 2759, 365, 521,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3533, 0, 6, 2759,
                                                                       203, 2786, 371, 539,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3587, 0, 6, 2786,
                                                                       212, 2813, 377, 557,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3641, 0, 6, 2813,
                                                                       221, 2840, 383, 575,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3695, 0, 6, 2840,
                                                                       230, 2867, 389, 593,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3749, 0, 6, 2867,
                                                                       239, 2894, 395, 611,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3803, 0, 6, 2921,
                                                                       275, 2948, 419, 665,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3857, 0, 6, 2948,
                                                                       284, 2975, 425, 683,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3911, 0, 6, 2975,
                                                                       293, 3002, 431, 701,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 3965, 0, 6, 3002,
                                                                       302, 3029, 437, 719,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 4019, 0, 6, 3029,
                                                                       311, 3056, 443, 737,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 4073, 0, 6, 3056,
                                                                       320, 3083, 449, 755,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 4127, 0, 6, 3083,
                                                                       329, 3110, 455, 773,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4181, 6, 359, 811,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4211, 6, 365, 821,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4241, 6, 371, 831,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4271, 6, 377, 841,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4301, 6, 383, 851,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4331, 6, 389, 861,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4361, 6, 395, 871,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4391, 6, 419, 901,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4421, 6, 425, 911,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4451, 6, 431, 921,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4481, 6, 437, 931,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4511, 6, 443, 941,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4541, 6, 449, 951,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 4571, 6, 455, 961,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4601, 0, 6, 3425,
                                                                       503, 3479, 811, 1031,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4691, 0, 6, 3479,
                                                                       521, 3533, 821, 1061,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4781, 0, 6, 3533,
                                                                       539, 3587, 831, 1091,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4871, 0, 6, 3587,
                                                                       557, 3641, 841, 1121,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 4961, 0, 6, 3641,
                                                                       575, 3695, 851, 1151,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 5051, 0, 6, 3695,
                                                                       593, 3749, 861, 1181,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 5141, 0, 6, 3803,
                                                                       665, 3857, 901, 1271,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 5231, 0, 6, 3857,
                                                                       683, 3911, 911, 1301,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 5321, 0, 6, 3911,
                                                                       701, 3965, 921, 1331,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 5411, 0, 6, 3965,
                                                                       719, 4019, 931, 1361,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 5501, 0, 6, 4019,
                                                                       737, 4073, 941, 1391,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 5591, 0, 6, 4073,
                                                                       755, 4127, 951, 1421,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5681, 6, 811,
                                                                       1481, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5726, 6, 821,
                                                                       1496, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5771, 6, 831,
                                                                       1511, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5816, 6, 841,
                                                                       1526, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5861, 6, 851,
                                                                       1541, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5906, 6, 861,
                                                                       1556, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5951, 6, 901,
                                                                       1601, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 5996, 6, 911,
                                                                       1616, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6041, 6, 921,
                                                                       1631, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6086, 6, 931,
                                                                       1646, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6131, 6, 941,
                                                                       1661, ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 6176, 6, 951,
                                                                       1676, ncols, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 6221, 0, 6, 4601,
                                                                       1031, 4691, 1481, 1781,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 6356, 0, 6, 4691,
                                                                       1061, 4781, 1496, 1826,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 6491, 0, 6, 4781,
                                                                       1091, 4871, 1511, 1871,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 6626, 0, 6, 4871,
                                                                       1121, 4961, 1526, 1916,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 6761, 0, 6, 4961,
                                                                       1151, 5051, 1541, 1961,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 6896, 0, 6, 5141,
                                                                       1271, 5231, 1601, 2096,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 7031, 0, 6, 5231,
                                                                       1301, 5321, 1616, 2141,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 7166, 0, 6, 5321,
                                                                       1331, 5411, 1631, 2186,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 7301, 0, 6, 5411,
                                                                       1361, 5501, 1646, 2231,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 7436, 0, 6, 5501,
                                                                       1391, 5591, 1661, 2276,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7571, 6, 10, 11,
                                                                       2321, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7577, 6, 11, 12,
                                                                       2324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7583, 6, 12, 13,
                                                                       2327, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7589, 6, 13, 14,
                                                                       2330, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7595, 6, 14, 15,
                                                                       2333, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7601, 6, 15, 16,
                                                                       2336, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7607, 6, 16, 17,
                                                                       2339, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7613, 6, 17, 18,
                                                                       2342, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7619, 6, 18, 19,
                                                                       2345, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7625, 6, 19, 20,
                                                                       2348, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7631, 6, 23, 24,
                                                                       2351, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7637, 6, 24, 25,
                                                                       2354, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7643, 6, 25, 26,
                                                                       2357, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7649, 6, 26, 27,
                                                                       2360, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7655, 6, 27, 28,
                                                                       2363, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7661, 6, 28, 29,
                                                                       2366, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7667, 6, 29, 30,
                                                                       2369, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7673, 6, 30, 31,
                                                                       2372, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7679, 6, 31, 32,
                                                                       2375, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 7685, 6, 32, 33,
                                                                       2378, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7691, 3, 6, 7571,
                                                                       2321, 7577, 2381, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7709, 3, 6, 7577,
                                                                       2324, 7583, 2390, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7727, 3, 6, 7583,
                                                                       2327, 7589, 2399, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7745, 3, 6, 7589,
                                                                       2330, 7595, 2408, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7763, 3, 6, 7595,
                                                                       2333, 7601, 2417, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7781, 3, 6, 7601,
                                                                       2336, 7607, 2426, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7799, 3, 6, 7607,
                                                                       2339, 7613, 2435, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7817, 3, 6, 7613,
                                                                       2342, 7619, 2444, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7835, 3, 6, 7619,
                                                                       2345, 7625, 2453, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7853, 3, 6, 7631,
                                                                       2351, 7637, 2462, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7871, 3, 6, 7637,
                                                                       2354, 7643, 2471, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7889, 3, 6, 7643,
                                                                       2357, 7649, 2480, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7907, 3, 6, 7649,
                                                                       2360, 7655, 2489, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7925, 3, 6, 7655,
                                                                       2363, 7661, 2498, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7943, 3, 6, 7661,
                                                                       2366, 7667, 2507, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7961, 3, 6, 7667,
                                                                       2369, 7673, 2516, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7979, 3, 6, 7673,
                                                                       2372, 7679, 2525, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 7997, 3, 6, 7679,
                                                                       2375, 7685, 2534, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8015, 0, 6, 7571,
                                                                       2321, 7577, 2543, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8033, 0, 6, 7577,
                                                                       2324, 7583, 2552, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8051, 0, 6, 7583,
                                                                       2327, 7589, 2561, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8069, 0, 6, 7589,
                                                                       2330, 7595, 2570, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8087, 0, 6, 7595,
                                                                       2333, 7601, 2579, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8105, 0, 6, 7601,
                                                                       2336, 7607, 2588, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8123, 0, 6, 7607,
                                                                       2339, 7613, 2597, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8141, 0, 6, 7613,
                                                                       2342, 7619, 2606, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8159, 0, 6, 7619,
                                                                       2345, 7625, 2615, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8177, 0, 6, 7631,
                                                                       2351, 7637, 2624, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8195, 0, 6, 7637,
                                                                       2354, 7643, 2633, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8213, 0, 6, 7643,
                                                                       2357, 7649, 2642, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8231, 0, 6, 7649,
                                                                       2360, 7655, 2651, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8249, 0, 6, 7655,
                                                                       2363, 7661, 2660, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8267, 0, 6, 7661,
                                                                       2366, 7667, 2669, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8285, 0, 6, 7667,
                                                                       2369, 7673, 2678, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8303, 0, 6, 7673,
                                                                       2372, 7679, 2687, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 8321, 0, 6, 7679,
                                                                       2375, 7685, 2696, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8339, 0, 3, 6,
                                                                       7691, 2381, 7709, 8015,
                                                                       2543, 8033, 167, 176,
                                                                       2705, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8393, 0, 3, 6,
                                                                       7709, 2390, 7727, 8033,
                                                                       2552, 8051, 176, 185,
                                                                       2732, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8447, 0, 3, 6,
                                                                       7727, 2399, 7745, 8051,
                                                                       2561, 8069, 185, 194,
                                                                       2759, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8501, 0, 3, 6,
                                                                       7745, 2408, 7763, 8069,
                                                                       2570, 8087, 194, 203,
                                                                       2786, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8555, 0, 3, 6,
                                                                       7763, 2417, 7781, 8087,
                                                                       2579, 8105, 203, 212,
                                                                       2813, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8609, 0, 3, 6,
                                                                       7781, 2426, 7799, 8105,
                                                                       2588, 8123, 212, 221,
                                                                       2840, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8663, 0, 3, 6,
                                                                       7799, 2435, 7817, 8123,
                                                                       2597, 8141, 221, 230,
                                                                       2867, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8717, 0, 3, 6,
                                                                       7817, 2444, 7835, 8141,
                                                                       2606, 8159, 230, 239,
                                                                       2894, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8771, 0, 3, 6,
                                                                       7853, 2462, 7871, 8177,
                                                                       2624, 8195, 257, 266,
                                                                       2921, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8825, 0, 3, 6,
                                                                       7871, 2471, 7889, 8195,
                                                                       2633, 8213, 266, 275,
                                                                       2948, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8879, 0, 3, 6,
                                                                       7889, 2480, 7907, 8213,
                                                                       2642, 8231, 275, 284,
                                                                       2975, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8933, 0, 3, 6,
                                                                       7907, 2489, 7925, 8231,
                                                                       2651, 8249, 284, 293,
                                                                       3002, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 8987, 0, 3, 6,
                                                                       7925, 2498, 7943, 8249,
                                                                       2660, 8267, 293, 302,
                                                                       3029, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9041, 0, 3, 6,
                                                                       7943, 2507, 7961, 8267,
                                                                       2669, 8285, 302, 311,
                                                                       3056, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9095, 0, 3, 6,
                                                                       7961, 2516, 7979, 8285,
                                                                       2678, 8303, 311, 320,
                                                                       3083, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 9149, 0, 3, 6,
                                                                       7979, 2525, 7997, 8303,
                                                                       2687, 8321, 320, 329,
                                                                       3110, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9203, 0, 6, 8015,
                                                                       2543, 8033, 347, 353,
                                                                       3137, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9239, 0, 6, 8033,
                                                                       2552, 8051, 353, 359,
                                                                       3155, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9275, 0, 6, 8051,
                                                                       2561, 8069, 359, 365,
                                                                       3173, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9311, 0, 6, 8069,
                                                                       2570, 8087, 365, 371,
                                                                       3191, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9347, 0, 6, 8087,
                                                                       2579, 8105, 371, 377,
                                                                       3209, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9383, 0, 6, 8105,
                                                                       2588, 8123, 377, 383,
                                                                       3227, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9419, 0, 6, 8123,
                                                                       2597, 8141, 383, 389,
                                                                       3245, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9455, 0, 6, 8141,
                                                                       2606, 8159, 389, 395,
                                                                       3263, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9491, 0, 6, 8177,
                                                                       2624, 8195, 407, 413,
                                                                       3281, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9527, 0, 6, 8195,
                                                                       2633, 8213, 413, 419,
                                                                       3299, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9563, 0, 6, 8213,
                                                                       2642, 8231, 419, 425,
                                                                       3317, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9599, 0, 6, 8231,
                                                                       2651, 8249, 425, 431,
                                                                       3335, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9635, 0, 6, 8249,
                                                                       2660, 8267, 431, 437,
                                                                       3353, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9671, 0, 6, 8267,
                                                                       2669, 8285, 437, 443,
                                                                       3371, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9707, 0, 6, 8285,
                                                                       2678, 8303, 443, 449,
                                                                       3389, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 9743, 0, 6, 8303,
                                                                       2687, 8321, 449, 455,
                                                                       3407, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 9779, 0, 3, 6,
                                                                       8339, 2705, 8393, 9203,
                                                                       3137, 9239, 467, 485,
                                                                       3425, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 9887, 0, 3, 6,
                                                                       8393, 2732, 8447, 9239,
                                                                       3155, 9275, 485, 503,
                                                                       3479, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 9995, 0, 3, 6,
                                                                       8447, 2759, 8501, 9275,
                                                                       3173, 9311, 503, 521,
                                                                       3533, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 10103, 0, 3, 6,
                                                                       8501, 2786, 8555, 9311,
                                                                       3191, 9347, 521, 539,
                                                                       3587, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 10211, 0, 3, 6,
                                                                       8555, 2813, 8609, 9347,
                                                                       3209, 9383, 539, 557,
                                                                       3641, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 10319, 0, 3, 6,
                                                                       8609, 2840, 8663, 9383,
                                                                       3227, 9419, 557, 575,
                                                                       3695, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 10427, 0, 3, 6,
                                                                       8663, 2867, 8717, 9419,
                                                                       3245, 9455, 575, 593,
                                                                       3749, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 10535, 0, 3, 6,
                                                                       8771, 2921, 8825, 9491,
                                                                       3281, 9527, 629, 647,
                                                                       3803, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 10643, 0, 3, 6,
                                                                       8825, 2948, 8879, 9527,
                                                                       3299, 9563, 647, 665,
                                                                       3857, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 10751, 0, 3, 6,
                                                                       8879, 2975, 8933, 9563,
                                                                       3317, 9599, 665, 683,
                                                                       3911, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 10859, 0, 3, 6,
                                                                       8933, 3002, 8987, 9599,
                                                                       3335, 9635, 683, 701,
                                                                       3965, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 10967, 0, 3, 6,
                                                                       8987, 3029, 9041, 9635,
                                                                       3353, 9671, 701, 719,
                                                                       4019, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 11075, 0, 3, 6,
                                                                       9041, 3056, 9095, 9671,
                                                                       3371, 9707, 719, 737,
                                                                       4073, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 11183, 0, 3, 6,
                                                                       9095, 3083, 9149, 9707,
                                                                       3389, 9743, 737, 755,
                                                                       4127, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11291, 0, 6, 9203,
                                                                       3137, 9239, 791, 801,
                                                                       4181, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11351, 0, 6, 9239,
                                                                       3155, 9275, 801, 811,
                                                                       4211, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11411, 0, 6, 9275,
                                                                       3173, 9311, 811, 821,
                                                                       4241, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11471, 0, 6, 9311,
                                                                       3191, 9347, 821, 831,
                                                                       4271, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11531, 0, 6, 9347,
                                                                       3209, 9383, 831, 841,
                                                                       4301, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11591, 0, 6, 9383,
                                                                       3227, 9419, 841, 851,
                                                                       4331, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11651, 0, 6, 9419,
                                                                       3245, 9455, 851, 861,
                                                                       4361, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11711, 0, 6, 9491,
                                                                       3281, 9527, 881, 891,
                                                                       4391, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11771, 0, 6, 9527,
                                                                       3299, 9563, 891, 901,
                                                                       4421, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11831, 0, 6, 9563,
                                                                       3317, 9599, 901, 911,
                                                                       4451, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11891, 0, 6, 9599,
                                                                       3335, 9635, 911, 921,
                                                                       4481, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 11951, 0, 6, 9635,
                                                                       3353, 9671, 921, 931,
                                                                       4511, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12011, 0, 6, 9671,
                                                                       3371, 9707, 931, 941,
                                                                       4541, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 12071, 0, 6, 9707,
                                                                       3389, 9743, 941, 951,
                                                                       4571, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 12131, 0, 3, 6,
                                                                       8339, 8393, 9779, 3425,
                                                                       9887, 11291, 4181, 11351,
                                                                       971, 1001, 4601, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 12311, 0, 3, 6,
                                                                       8393, 8447, 9887, 3479,
                                                                       9995, 11351, 4211, 11411,
                                                                       1001, 1031, 4691, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 12491, 0, 3, 6,
                                                                       8447, 8501, 9995, 3533,
                                                                       10103, 11411, 4241, 11471,
                                                                       1031, 1061, 4781, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 12671, 0, 3, 6,
                                                                       8501, 8555, 10103, 3587,
                                                                       10211, 11471, 4271, 11531,
                                                                       1061, 1091, 4871, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 12851, 0, 3, 6,
                                                                       8555, 8609, 10211, 3641,
                                                                       10319, 11531, 4301, 11591,
                                                                       1091, 1121, 4961, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 13031, 0, 3, 6,
                                                                       8609, 8663, 10319, 3695,
                                                                       10427, 11591, 4331, 11651,
                                                                       1121, 1151, 5051, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 13211, 0, 3, 6,
                                                                       8771, 8825, 10535, 3803,
                                                                       10643, 11711, 4391, 11771,
                                                                       1211, 1241, 5141, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 13391, 0, 3, 6,
                                                                       8825, 8879, 10643, 3857,
                                                                       10751, 11771, 4421, 11831,
                                                                       1241, 1271, 5231, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 13571, 0, 3, 6,
                                                                       8879, 8933, 10751, 3911,
                                                                       10859, 11831, 4451, 11891,
                                                                       1271, 1301, 5321, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 13751, 0, 3, 6,
                                                                       8933, 8987, 10859, 3965,
                                                                       10967, 11891, 4481, 11951,
                                                                       1301, 1331, 5411, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 13931, 0, 3, 6,
                                                                       8987, 9041, 10967, 4019,
                                                                       11075, 11951, 4511, 12011,
                                                                       1331, 1361, 5501, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 14111, 0, 3, 6,
                                                                       9041, 9095, 11075, 4073,
                                                                       11183, 12011, 4541, 12071,
                                                                       1361, 1391, 5591, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14291, 0, 6,
                                                                       11291, 4181, 11351, 1451,
                                                                       1466, 5681, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14381, 0, 6,
                                                                       11351, 4211, 11411, 1466,
                                                                       1481, 5726, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14471, 0, 6,
                                                                       11411, 4241, 11471, 1481,
                                                                       1496, 5771, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14561, 0, 6,
                                                                       11471, 4271, 11531, 1496,
                                                                       1511, 5816, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14651, 0, 6,
                                                                       11531, 4301, 11591, 1511,
                                                                       1526, 5861, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14741, 0, 6,
                                                                       11591, 4331, 11651, 1526,
                                                                       1541, 5906, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14831, 0, 6,
                                                                       11711, 4391, 11771, 1571,
                                                                       1586, 5951, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 14921, 0, 6,
                                                                       11771, 4421, 11831, 1586,
                                                                       1601, 5996, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15011, 0, 6,
                                                                       11831, 4451, 11891, 1601,
                                                                       1616, 6041, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15101, 0, 6,
                                                                       11891, 4481, 11951, 1616,
                                                                       1631, 6086, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15191, 0, 6,
                                                                       11951, 4511, 12011, 1631,
                                                                       1646, 6131, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 15281, 0, 6,
                                                                       12011, 4541, 12071, 1646,
                                                                       1661, 6176, ncols, gamma,
                                                                       p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 15371, 0, 3, 6,
                                                                       9779, 9887, 12131, 4601,
                                                                       12311, 14291, 5681, 14381,
                                                                       1691, 1736, 6221, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 15641, 0, 3, 6,
                                                                       9887, 9995, 12311, 4691,
                                                                       12491, 14381, 5726, 14471,
                                                                       1736, 1781, 6356, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 15911, 0, 3, 6,
                                                                       9995, 10103, 12491, 4781,
                                                                       12671, 14471, 5771, 14561,
                                                                       1781, 1826, 6491, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 16181, 0, 3, 6,
                                                                       10103, 10211, 12671, 4871,
                                                                       12851, 14561, 5816, 14651,
                                                                       1826, 1871, 6626, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 16451, 0, 3, 6,
                                                                       10211, 10319, 12851, 4961,
                                                                       13031, 14651, 5861, 14741,
                                                                       1871, 1916, 6761, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 16721, 0, 3, 6,
                                                                       10535, 10643, 13211, 5141,
                                                                       13391, 14831, 5951, 14921,
                                                                       2006, 2051, 6896, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 16991, 0, 3, 6,
                                                                       10643, 10751, 13391, 5231,
                                                                       13571, 14921, 5996, 15011,
                                                                       2051, 2096, 7031, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 17261, 0, 3, 6,
                                                                       10751, 10859, 13571, 5321,
                                                                       13751, 15011, 6041, 15101,
                                                                       2096, 2141, 7166, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 17531, 0, 3, 6,
                                                                       10859, 10967, 13751, 5411,
                                                                       13931, 15101, 6086, 15191,
                                                                       2141, 2186, 7301, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 17801, 0, 3, 6,
                                                                       10967, 11075, 13931, 5501,
                                                                       14111, 15191, 6131, 15281,
                                                                       2186, 2231, 7436, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18071, 6, 2321,
                                                                       2324, 7583, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18081, 6, 2324,
                                                                       2327, 7589, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18091, 6, 2327,
                                                                       2330, 7595, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18101, 6, 2330,
                                                                       2333, 7601, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18111, 6, 2333,
                                                                       2336, 7607, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18121, 6, 2336,
                                                                       2339, 7613, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18131, 6, 2339,
                                                                       2342, 7619, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18141, 6, 2342,
                                                                       2345, 7625, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18151, 6, 2351,
                                                                       2354, 7643, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18161, 6, 2354,
                                                                       2357, 7649, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18171, 6, 2357,
                                                                       2360, 7655, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18181, 6, 2360,
                                                                       2363, 7661, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18191, 6, 2363,
                                                                       2366, 7667, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18201, 6, 2366,
                                                                       2369, 7673, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18211, 6, 2369,
                                                                       2372, 7679, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 18221, 6, 2372,
                                                                       2375, 7685, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18231, 3, 6,
                                                                       18071, 7583, 18081, 7727,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18261, 3, 6,
                                                                       18081, 7589, 18091, 7745,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18291, 3, 6,
                                                                       18091, 7595, 18101, 7763,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18321, 3, 6,
                                                                       18101, 7601, 18111, 7781,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18351, 3, 6,
                                                                       18111, 7607, 18121, 7799,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18381, 3, 6,
                                                                       18121, 7613, 18131, 7817,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18411, 3, 6,
                                                                       18131, 7619, 18141, 7835,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18441, 3, 6,
                                                                       18151, 7643, 18161, 7889,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18471, 3, 6,
                                                                       18161, 7649, 18171, 7907,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18501, 3, 6,
                                                                       18171, 7655, 18181, 7925,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18531, 3, 6,
                                                                       18181, 7661, 18191, 7943,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18561, 3, 6,
                                                                       18191, 7667, 18201, 7961,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18591, 3, 6,
                                                                       18201, 7673, 18211, 7979,
                                                                       ncols, gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 18621, 3, 6,
                                                                       18211, 7679, 18221, 7997,
                                                                       ncols, gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18651, 0, 6,
                                                                       18071, 7583, 18081, 2543,
                                                                       2552, 8051, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18681, 0, 6,
                                                                       18081, 7589, 18091, 2552,
                                                                       2561, 8069, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18711, 0, 6,
                                                                       18091, 7595, 18101, 2561,
                                                                       2570, 8087, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18741, 0, 6,
                                                                       18101, 7601, 18111, 2570,
                                                                       2579, 8105, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18771, 0, 6,
                                                                       18111, 7607, 18121, 2579,
                                                                       2588, 8123, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18801, 0, 6,
                                                                       18121, 7613, 18131, 2588,
                                                                       2597, 8141, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18831, 0, 6,
                                                                       18131, 7619, 18141, 2597,
                                                                       2606, 8159, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18861, 0, 6,
                                                                       18151, 7643, 18161, 2624,
                                                                       2633, 8213, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18891, 0, 6,
                                                                       18161, 7649, 18171, 2633,
                                                                       2642, 8231, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18921, 0, 6,
                                                                       18171, 7655, 18181, 2642,
                                                                       2651, 8249, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18951, 0, 6,
                                                                       18181, 7661, 18191, 2651,
                                                                       2660, 8267, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 18981, 0, 6,
                                                                       18191, 7667, 18201, 2660,
                                                                       2669, 8285, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19011, 0, 6,
                                                                       18201, 7673, 18211, 2669,
                                                                       2678, 8303, ncols, gamma,
                                                                       p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 19041, 0, 6,
                                                                       18211, 7679, 18221, 2678,
                                                                       2687, 8321, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 19071, 0, 3, 6,
                                                                       18231, 7727, 18261, 18651,
                                                                       8051, 18681, 2705, 2732,
                                                                       8447, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 19161, 0, 3, 6,
                                                                       18261, 7745, 18291, 18681,
                                                                       8069, 18711, 2732, 2759,
                                                                       8501, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 19251, 0, 3, 6,
                                                                       18291, 7763, 18321, 18711,
                                                                       8087, 18741, 2759, 2786,
                                                                       8555, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 19341, 0, 3, 6,
                                                                       18321, 7781, 18351, 18741,
                                                                       8105, 18771, 2786, 2813,
                                                                       8609, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 19431, 0, 3, 6,
                                                                       18351, 7799, 18381, 18771,
                                                                       8123, 18801, 2813, 2840,
                                                                       8663, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 19521, 0, 3, 6,
                                                                       18381, 7817, 18411, 18801,
                                                                       8141, 18831, 2840, 2867,
                                                                       8717, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 19611, 0, 3, 6,
                                                                       18441, 7889, 18471, 18861,
                                                                       8213, 18891, 2921, 2948,
                                                                       8879, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 19701, 0, 3, 6,
                                                                       18471, 7907, 18501, 18891,
                                                                       8231, 18921, 2948, 2975,
                                                                       8933, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 19791, 0, 3, 6,
                                                                       18501, 7925, 18531, 18921,
                                                                       8249, 18951, 2975, 3002,
                                                                       8987, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 19881, 0, 3, 6,
                                                                       18531, 7943, 18561, 18951,
                                                                       8267, 18981, 3002, 3029,
                                                                       9041, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 19971, 0, 3, 6,
                                                                       18561, 7961, 18591, 18981,
                                                                       8285, 19011, 3029, 3056,
                                                                       9095, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 20061, 0, 3, 6,
                                                                       18591, 7979, 18621, 19011,
                                                                       8303, 19041, 3056, 3083,
                                                                       9149, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20151, 0, 6,
                                                                       18651, 8051, 18681, 3137,
                                                                       3155, 9275, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20211, 0, 6,
                                                                       18681, 8069, 18711, 3155,
                                                                       3173, 9311, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20271, 0, 6,
                                                                       18711, 8087, 18741, 3173,
                                                                       3191, 9347, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20331, 0, 6,
                                                                       18741, 8105, 18771, 3191,
                                                                       3209, 9383, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20391, 0, 6,
                                                                       18771, 8123, 18801, 3209,
                                                                       3227, 9419, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20451, 0, 6,
                                                                       18801, 8141, 18831, 3227,
                                                                       3245, 9455, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20511, 0, 6,
                                                                       18861, 8213, 18891, 3281,
                                                                       3299, 9563, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20571, 0, 6,
                                                                       18891, 8231, 18921, 3299,
                                                                       3317, 9599, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20631, 0, 6,
                                                                       18921, 8249, 18951, 3317,
                                                                       3335, 9635, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20691, 0, 6,
                                                                       18951, 8267, 18981, 3335,
                                                                       3353, 9671, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20751, 0, 6,
                                                                       18981, 8285, 19011, 3353,
                                                                       3371, 9707, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 20811, 0, 6,
                                                                       19011, 8303, 19041, 3371,
                                                                       3389, 9743, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 20871, 0, 3, 6,
                                                                       19071, 8447, 19161, 20151,
                                                                       9275, 20211, 3425, 3479,
                                                                       9995, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 21051, 0, 3, 6,
                                                                       19161, 8501, 19251, 20211,
                                                                       9311, 20271, 3479, 3533,
                                                                       10103, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 21231, 0, 3, 6,
                                                                       19251, 8555, 19341, 20271,
                                                                       9347, 20331, 3533, 3587,
                                                                       10211, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 21411, 0, 3, 6,
                                                                       19341, 8609, 19431, 20331,
                                                                       9383, 20391, 3587, 3641,
                                                                       10319, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 21591, 0, 3, 6,
                                                                       19431, 8663, 19521, 20391,
                                                                       9419, 20451, 3641, 3695,
                                                                       10427, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 21771, 0, 3, 6,
                                                                       19611, 8879, 19701, 20511,
                                                                       9563, 20571, 3803, 3857,
                                                                       10751, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 21951, 0, 3, 6,
                                                                       19701, 8933, 19791, 20571,
                                                                       9599, 20631, 3857, 3911,
                                                                       10859, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 22131, 0, 3, 6,
                                                                       19791, 8987, 19881, 20631,
                                                                       9635, 20691, 3911, 3965,
                                                                       10967, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 22311, 0, 3, 6,
                                                                       19881, 9041, 19971, 20691,
                                                                       9671, 20751, 3965, 4019,
                                                                       11075, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 22491, 0, 3, 6,
                                                                       19971, 9095, 20061, 20751,
                                                                       9707, 20811, 4019, 4073,
                                                                       11183, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22671, 0, 6,
                                                                       20151, 9275, 20211, 4181,
                                                                       4211, 11411, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22771, 0, 6,
                                                                       20211, 9311, 20271, 4211,
                                                                       4241, 11471, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22871, 0, 6,
                                                                       20271, 9347, 20331, 4241,
                                                                       4271, 11531, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 22971, 0, 6,
                                                                       20331, 9383, 20391, 4271,
                                                                       4301, 11591, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23071, 0, 6,
                                                                       20391, 9419, 20451, 4301,
                                                                       4331, 11651, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23171, 0, 6,
                                                                       20511, 9563, 20571, 4391,
                                                                       4421, 11831, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23271, 0, 6,
                                                                       20571, 9599, 20631, 4421,
                                                                       4451, 11891, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23371, 0, 6,
                                                                       20631, 9635, 20691, 4451,
                                                                       4481, 11951, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23471, 0, 6,
                                                                       20691, 9671, 20751, 4481,
                                                                       4511, 12011, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 23571, 0, 6,
                                                                       20751, 9707, 20811, 4511,
                                                                       4541, 12071, ncols, gamma,
                                                                       p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 23671, 0, 3, 6,
                                                                       19071, 19161, 20871, 9995,
                                                                       21051, 22671, 11411,
                                                                       22771, 4601, 4691, 12491,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 23971, 0, 3, 6,
                                                                       19161, 19251, 21051,
                                                                       10103, 21231, 22771,
                                                                       11471, 22871, 4691, 4781,
                                                                       12671, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 24271, 0, 3, 6,
                                                                       19251, 19341, 21231,
                                                                       10211, 21411, 22871,
                                                                       11531, 22971, 4781, 4871,
                                                                       12851, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 24571, 0, 3, 6,
                                                                       19341, 19431, 21411,
                                                                       10319, 21591, 22971,
                                                                       11591, 23071, 4871, 4961,
                                                                       13031, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 24871, 0, 3, 6,
                                                                       19611, 19701, 21771,
                                                                       10751, 21951, 23171,
                                                                       11831, 23271, 5141, 5231,
                                                                       13571, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 25171, 0, 3, 6,
                                                                       19701, 19791, 21951,
                                                                       10859, 22131, 23271,
                                                                       11891, 23371, 5231, 5321,
                                                                       13751, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 25471, 0, 3, 6,
                                                                       19791, 19881, 22131,
                                                                       10967, 22311, 23371,
                                                                       11951, 23471, 5321, 5411,
                                                                       13931, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 25771, 0, 3, 6,
                                                                       19881, 19971, 22311,
                                                                       11075, 22491, 23471,
                                                                       12011, 23571, 5411, 5501,
                                                                       14111, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26071, 0, 6,
                                                                       22671, 11411, 22771, 5681,
                                                                       5726, 14471, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26221, 0, 6,
                                                                       22771, 11471, 22871, 5726,
                                                                       5771, 14561, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26371, 0, 6,
                                                                       22871, 11531, 22971, 5771,
                                                                       5816, 14651, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26521, 0, 6,
                                                                       22971, 11591, 23071, 5816,
                                                                       5861, 14741, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26671, 0, 6,
                                                                       23171, 11831, 23271, 5951,
                                                                       5996, 15011, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26821, 0, 6,
                                                                       23271, 11891, 23371, 5996,
                                                                       6041, 15101, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 26971, 0, 6,
                                                                       23371, 11951, 23471, 6041,
                                                                       6086, 15191, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 27121, 0, 6,
                                                                       23471, 12011, 23571, 6086,
                                                                       6131, 15281, ncols, gamma,
                                                                       p, q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 27271, 0, 3, 6,
                                                                       20871, 21051, 23671,
                                                                       12491, 23971, 26071,
                                                                       14471, 26221, 6221, 6356,
                                                                       15911, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 27721, 0, 3, 6,
                                                                       21051, 21231, 23971,
                                                                       12671, 24271, 26221,
                                                                       14561, 26371, 6356, 6491,
                                                                       16181, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 28171, 0, 3, 6,
                                                                       21231, 21411, 24271,
                                                                       12851, 24571, 26371,
                                                                       14651, 26521, 6491, 6626,
                                                                       16451, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 28621, 0, 3, 6,
                                                                       21771, 21951, 24871,
                                                                       13571, 25171, 26671,
                                                                       15011, 26821, 6896, 7031,
                                                                       17261, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 29071, 0, 3, 6,
                                                                       21951, 22131, 25171,
                                                                       13751, 25471, 26821,
                                                                       15101, 26971, 7031, 7166,
                                                                       17531, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 29521, 0, 3, 6,
                                                                       22131, 22311, 25471,
                                                                       13931, 25771, 26971,
                                                                       15191, 27121, 7166, 7301,
                                                                       17801, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 29971, 6, 7571,
                                                                       7577, 18071, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 29986, 6, 7577,
                                                                       7583, 18081, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30001, 6, 7583,
                                                                       7589, 18091, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30016, 6, 7589,
                                                                       7595, 18101, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30031, 6, 7595,
                                                                       7601, 18111, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30046, 6, 7601,
                                                                       7607, 18121, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30061, 6, 7607,
                                                                       7613, 18131, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30076, 6, 7613,
                                                                       7619, 18141, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30091, 6, 7631,
                                                                       7637, 18151, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30106, 6, 7637,
                                                                       7643, 18161, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30121, 6, 7643,
                                                                       7649, 18171, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30136, 6, 7649,
                                                                       7655, 18181, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30151, 6, 7655,
                                                                       7661, 18191, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30166, 6, 7661,
                                                                       7667, 18201, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30181, 6, 7667,
                                                                       7673, 18211, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 30196, 6, 7673,
                                                                       7679, 18221, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30211, 3, 6,
                                                                       29971, 18071, 29986, 7691,
                                                                       7709, 18231, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30256, 3, 6,
                                                                       29986, 18081, 30001, 7709,
                                                                       7727, 18261, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30301, 3, 6,
                                                                       30001, 18091, 30016, 7727,
                                                                       7745, 18291, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30346, 3, 6,
                                                                       30016, 18101, 30031, 7745,
                                                                       7763, 18321, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30391, 3, 6,
                                                                       30031, 18111, 30046, 7763,
                                                                       7781, 18351, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30436, 3, 6,
                                                                       30046, 18121, 30061, 7781,
                                                                       7799, 18381, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30481, 3, 6,
                                                                       30061, 18131, 30076, 7799,
                                                                       7817, 18411, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30526, 3, 6,
                                                                       30091, 18151, 30106, 7853,
                                                                       7871, 18441, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30571, 3, 6,
                                                                       30106, 18161, 30121, 7871,
                                                                       7889, 18471, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30616, 3, 6,
                                                                       30121, 18171, 30136, 7889,
                                                                       7907, 18501, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30661, 3, 6,
                                                                       30136, 18181, 30151, 7907,
                                                                       7925, 18531, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30706, 3, 6,
                                                                       30151, 18191, 30166, 7925,
                                                                       7943, 18561, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30751, 3, 6,
                                                                       30166, 18201, 30181, 7943,
                                                                       7961, 18591, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 30796, 3, 6,
                                                                       30181, 18211, 30196, 7961,
                                                                       7979, 18621, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 30841, 0, 6,
                                                                       29971, 18071, 29986, 8015,
                                                                       8033, 18651, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 30886, 0, 6,
                                                                       29986, 18081, 30001, 8033,
                                                                       8051, 18681, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 30931, 0, 6,
                                                                       30001, 18091, 30016, 8051,
                                                                       8069, 18711, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 30976, 0, 6,
                                                                       30016, 18101, 30031, 8069,
                                                                       8087, 18741, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 31021, 0, 6,
                                                                       30031, 18111, 30046, 8087,
                                                                       8105, 18771, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 31066, 0, 6,
                                                                       30046, 18121, 30061, 8105,
                                                                       8123, 18801, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 31111, 0, 6,
                                                                       30061, 18131, 30076, 8123,
                                                                       8141, 18831, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 31156, 0, 6,
                                                                       30091, 18151, 30106, 8177,
                                                                       8195, 18861, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 31201, 0, 6,
                                                                       30106, 18161, 30121, 8195,
                                                                       8213, 18891, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 31246, 0, 6,
                                                                       30121, 18171, 30136, 8213,
                                                                       8231, 18921, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 31291, 0, 6,
                                                                       30136, 18181, 30151, 8231,
                                                                       8249, 18951, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 31336, 0, 6,
                                                                       30151, 18191, 30166, 8249,
                                                                       8267, 18981, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 31381, 0, 6,
                                                                       30166, 18201, 30181, 8267,
                                                                       8285, 19011, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 31426, 0, 6,
                                                                       30181, 18211, 30196, 8285,
                                                                       8303, 19041, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 31471, 0, 3, 6,
                                                                       30211, 18231, 30256,
                                                                       30841, 18651, 30886, 8339,
                                                                       8393, 19071, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 31606, 0, 3, 6,
                                                                       30256, 18261, 30301,
                                                                       30886, 18681, 30931, 8393,
                                                                       8447, 19161, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 31741, 0, 3, 6,
                                                                       30301, 18291, 30346,
                                                                       30931, 18711, 30976, 8447,
                                                                       8501, 19251, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 31876, 0, 3, 6,
                                                                       30346, 18321, 30391,
                                                                       30976, 18741, 31021, 8501,
                                                                       8555, 19341, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 32011, 0, 3, 6,
                                                                       30391, 18351, 30436,
                                                                       31021, 18771, 31066, 8555,
                                                                       8609, 19431, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 32146, 0, 3, 6,
                                                                       30436, 18381, 30481,
                                                                       31066, 18801, 31111, 8609,
                                                                       8663, 19521, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 32281, 0, 3, 6,
                                                                       30526, 18441, 30571,
                                                                       31156, 18861, 31201, 8771,
                                                                       8825, 19611, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 32416, 0, 3, 6,
                                                                       30571, 18471, 30616,
                                                                       31201, 18891, 31246, 8825,
                                                                       8879, 19701, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 32551, 0, 3, 6,
                                                                       30616, 18501, 30661,
                                                                       31246, 18921, 31291, 8879,
                                                                       8933, 19791, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 32686, 0, 3, 6,
                                                                       30661, 18531, 30706,
                                                                       31291, 18951, 31336, 8933,
                                                                       8987, 19881, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 32821, 0, 3, 6,
                                                                       30706, 18561, 30751,
                                                                       31336, 18981, 31381, 8987,
                                                                       9041, 19971, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 32956, 0, 3, 6,
                                                                       30751, 18591, 30796,
                                                                       31381, 19011, 31426, 9041,
                                                                       9095, 20061, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 33091, 0, 6,
                                                                       30841, 18651, 30886, 9203,
                                                                       9239, 20151, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 33181, 0, 6,
                                                                       30886, 18681, 30931, 9239,
                                                                       9275, 20211, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 33271, 0, 6,
                                                                       30931, 18711, 30976, 9275,
                                                                       9311, 20271, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 33361, 0, 6,
                                                                       30976, 18741, 31021, 9311,
                                                                       9347, 20331, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 33451, 0, 6,
                                                                       31021, 18771, 31066, 9347,
                                                                       9383, 20391, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 33541, 0, 6,
                                                                       31066, 18801, 31111, 9383,
                                                                       9419, 20451, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 33631, 0, 6,
                                                                       31156, 18861, 31201, 9491,
                                                                       9527, 20511, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 33721, 0, 6,
                                                                       31201, 18891, 31246, 9527,
                                                                       9563, 20571, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 33811, 0, 6,
                                                                       31246, 18921, 31291, 9563,
                                                                       9599, 20631, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 33901, 0, 6,
                                                                       31291, 18951, 31336, 9599,
                                                                       9635, 20691, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 33991, 0, 6,
                                                                       31336, 18981, 31381, 9635,
                                                                       9671, 20751, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 34081, 0, 6,
                                                                       31381, 19011, 31426, 9671,
                                                                       9707, 20811, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 34171, 0, 3, 6,
                                                                       31471, 19071, 31606,
                                                                       33091, 20151, 33181, 9779,
                                                                       9887, 20871, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 34441, 0, 3, 6,
                                                                       31606, 19161, 31741,
                                                                       33181, 20211, 33271, 9887,
                                                                       9995, 21051, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 34711, 0, 3, 6,
                                                                       31741, 19251, 31876,
                                                                       33271, 20271, 33361, 9995,
                                                                       10103, 21231, ncols,
                                                                       gamma, p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 34981, 0, 3, 6,
                                                                       31876, 19341, 32011,
                                                                       33361, 20331, 33451,
                                                                       10103, 10211, 21411,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 35251, 0, 3, 6,
                                                                       32011, 19431, 32146,
                                                                       33451, 20391, 33541,
                                                                       10211, 10319, 21591,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 35521, 0, 3, 6,
                                                                       32281, 19611, 32416,
                                                                       33631, 20511, 33721,
                                                                       10535, 10643, 21771,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 35791, 0, 3, 6,
                                                                       32416, 19701, 32551,
                                                                       33721, 20571, 33811,
                                                                       10643, 10751, 21951,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 36061, 0, 3, 6,
                                                                       32551, 19791, 32686,
                                                                       33811, 20631, 33901,
                                                                       10751, 10859, 22131,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 36331, 0, 3, 6,
                                                                       32686, 19881, 32821,
                                                                       33901, 20691, 33991,
                                                                       10859, 10967, 22311,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 36601, 0, 3, 6,
                                                                       32821, 19971, 32956,
                                                                       33991, 20751, 34081,
                                                                       10967, 11075, 22491,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 36871, 0, 6,
                                                                       33091, 20151, 33181,
                                                                       11291, 11351, 22671,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37021, 0, 6,
                                                                       33181, 20211, 33271,
                                                                       11351, 11411, 22771,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37171, 0, 6,
                                                                       33271, 20271, 33361,
                                                                       11411, 11471, 22871,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37321, 0, 6,
                                                                       33361, 20331, 33451,
                                                                       11471, 11531, 22971,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37471, 0, 6,
                                                                       33451, 20391, 33541,
                                                                       11531, 11591, 23071,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37621, 0, 6,
                                                                       33631, 20511, 33721,
                                                                       11711, 11771, 23171,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37771, 0, 6,
                                                                       33721, 20571, 33811,
                                                                       11771, 11831, 23271,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 37921, 0, 6,
                                                                       33811, 20631, 33901,
                                                                       11831, 11891, 23371,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38071, 0, 6,
                                                                       33901, 20691, 33991,
                                                                       11891, 11951, 23471,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 38221, 0, 6,
                                                                       33991, 20751, 34081,
                                                                       11951, 12011, 23571,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 38371, 0, 3, 6,
                                                                       31471, 31606, 34171,
                                                                       20871, 34441, 36871,
                                                                       22671, 37021, 12131,
                                                                       12311, 23671, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 38821, 0, 3, 6,
                                                                       31606, 31741, 34441,
                                                                       21051, 34711, 37021,
                                                                       22771, 37171, 12311,
                                                                       12491, 23971, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 39271, 0, 3, 6,
                                                                       31741, 31876, 34711,
                                                                       21231, 34981, 37171,
                                                                       22871, 37321, 12491,
                                                                       12671, 24271, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 39721, 0, 3, 6,
                                                                       31876, 32011, 34981,
                                                                       21411, 35251, 37321,
                                                                       22971, 37471, 12671,
                                                                       12851, 24571, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 40171, 0, 3, 6,
                                                                       32281, 32416, 35521,
                                                                       21771, 35791, 37621,
                                                                       23171, 37771, 13211,
                                                                       13391, 24871, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 40621, 0, 3, 6,
                                                                       32416, 32551, 35791,
                                                                       21951, 36061, 37771,
                                                                       23271, 37921, 13391,
                                                                       13571, 25171, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 41071, 0, 3, 6,
                                                                       32551, 32686, 36061,
                                                                       22131, 36331, 37921,
                                                                       23371, 38071, 13571,
                                                                       13751, 25471, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 41521, 0, 3, 6,
                                                                       32686, 32821, 36331,
                                                                       22311, 36601, 38071,
                                                                       23471, 38221, 13751,
                                                                       13931, 25771, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 41971, 0, 6,
                                                                       36871, 22671, 37021,
                                                                       14291, 14381, 26071,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 42196, 0, 6,
                                                                       37021, 22771, 37171,
                                                                       14381, 14471, 26221,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 42421, 0, 6,
                                                                       37171, 22871, 37321,
                                                                       14471, 14561, 26371,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 42646, 0, 6,
                                                                       37321, 22971, 37471,
                                                                       14561, 14651, 26521,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 42871, 0, 6,
                                                                       37621, 23171, 37771,
                                                                       14831, 14921, 26671,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 43096, 0, 6,
                                                                       37771, 23271, 37921,
                                                                       14921, 15011, 26821,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 43321, 0, 6,
                                                                       37921, 23371, 38071,
                                                                       15011, 15101, 26971,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 43546, 0, 6,
                                                                       38071, 23471, 38221,
                                                                       15101, 15191, 27121,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 43771, 0, 3, 6,
                                                                       34171, 34441, 38371,
                                                                       23671, 38821, 41971,
                                                                       26071, 42196, 15371,
                                                                       15641, 27271, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 44446, 0, 3, 6,
                                                                       34441, 34711, 38821,
                                                                       23971, 39271, 42196,
                                                                       26221, 42421, 15641,
                                                                       15911, 27721, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 45121, 0, 3, 6,
                                                                       34711, 34981, 39271,
                                                                       24271, 39721, 42421,
                                                                       26371, 42646, 15911,
                                                                       16181, 28171, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 45796, 0, 3, 6,
                                                                       35521, 35791, 40171,
                                                                       24871, 40621, 42871,
                                                                       26671, 43096, 16721,
                                                                       16991, 28621, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 46471, 0, 3, 6,
                                                                       35791, 36061, 40621,
                                                                       25171, 41071, 43096,
                                                                       26821, 43321, 16991,
                                                                       17261, 29071, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 47146, 0, 3, 6,
                                                                       36061, 36331, 41071,
                                                                       25471, 41521, 43321,
                                                                       26971, 43546, 17261,
                                                                       17531, 29521, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47821, 6, 18071,
                                                                       18081, 30001, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47842, 6, 18081,
                                                                       18091, 30016, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47863, 6, 18091,
                                                                       18101, 30031, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47884, 6, 18101,
                                                                       18111, 30046, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47905, 6, 18111,
                                                                       18121, 30061, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47926, 6, 18121,
                                                                       18131, 30076, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47947, 6, 18151,
                                                                       18161, 30121, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47968, 6, 18161,
                                                                       18171, 30136, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 47989, 6, 18171,
                                                                       18181, 30151, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48010, 6, 18181,
                                                                       18191, 30166, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48031, 6, 18191,
                                                                       18201, 30181, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 48052, 6, 18201,
                                                                       18211, 30196, ncols,
                                                                       gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48073, 3, 6,
                                                                       47821, 30001, 47842,
                                                                       18231, 18261, 30301,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48136, 3, 6,
                                                                       47842, 30016, 47863,
                                                                       18261, 18291, 30346,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48199, 3, 6,
                                                                       47863, 30031, 47884,
                                                                       18291, 18321, 30391,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48262, 3, 6,
                                                                       47884, 30046, 47905,
                                                                       18321, 18351, 30436,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48325, 3, 6,
                                                                       47905, 30061, 47926,
                                                                       18351, 18381, 30481,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48388, 3, 6,
                                                                       47947, 30121, 47968,
                                                                       18441, 18471, 30616,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48451, 3, 6,
                                                                       47968, 30136, 47989,
                                                                       18471, 18501, 30661,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48514, 3, 6,
                                                                       47989, 30151, 48010,
                                                                       18501, 18531, 30706,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48577, 3, 6,
                                                                       48010, 30166, 48031,
                                                                       18531, 18561, 30751,
                                                                       ncols, gamma, p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 48640, 3, 6,
                                                                       48031, 30181, 48052,
                                                                       18561, 18591, 30796,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 48703, 0, 6,
                                                                       47821, 30001, 47842,
                                                                       18651, 18681, 30931,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 48766, 0, 6,
                                                                       47842, 30016, 47863,
                                                                       18681, 18711, 30976,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 48829, 0, 6,
                                                                       47863, 30031, 47884,
                                                                       18711, 18741, 31021,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 48892, 0, 6,
                                                                       47884, 30046, 47905,
                                                                       18741, 18771, 31066,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 48955, 0, 6,
                                                                       47905, 30061, 47926,
                                                                       18771, 18801, 31111,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49018, 0, 6,
                                                                       47947, 30121, 47968,
                                                                       18861, 18891, 31246,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49081, 0, 6,
                                                                       47968, 30136, 47989,
                                                                       18891, 18921, 31291,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49144, 0, 6,
                                                                       47989, 30151, 48010,
                                                                       18921, 18951, 31336,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49207, 0, 6,
                                                                       48010, 30166, 48031,
                                                                       18951, 18981, 31381,
                                                                       ncols, gamma, p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 49270, 0, 6,
                                                                       48031, 30181, 48052,
                                                                       18981, 19011, 31426,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 49333, 0, 3, 6,
                                                                       48073, 30301, 48136,
                                                                       48703, 30931, 48766,
                                                                       19071, 19161, 31741,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 49522, 0, 3, 6,
                                                                       48136, 30346, 48199,
                                                                       48766, 30976, 48829,
                                                                       19161, 19251, 31876,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 49711, 0, 3, 6,
                                                                       48199, 30391, 48262,
                                                                       48829, 31021, 48892,
                                                                       19251, 19341, 32011,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 49900, 0, 3, 6,
                                                                       48262, 30436, 48325,
                                                                       48892, 31066, 48955,
                                                                       19341, 19431, 32146,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 50089, 0, 3, 6,
                                                                       48388, 30616, 48451,
                                                                       49018, 31246, 49081,
                                                                       19611, 19701, 32551,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 50278, 0, 3, 6,
                                                                       48451, 30661, 48514,
                                                                       49081, 31291, 49144,
                                                                       19701, 19791, 32686,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 50467, 0, 3, 6,
                                                                       48514, 30706, 48577,
                                                                       49144, 31336, 49207,
                                                                       19791, 19881, 32821,
                                                                       ncols, gamma, p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 50656, 0, 3, 6,
                                                                       48577, 30751, 48640,
                                                                       49207, 31381, 49270,
                                                                       19881, 19971, 32956,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 50845, 0, 6,
                                                                       48703, 30931, 48766,
                                                                       20151, 20211, 33271,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 50971, 0, 6,
                                                                       48766, 30976, 48829,
                                                                       20211, 20271, 33361,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 51097, 0, 6,
                                                                       48829, 31021, 48892,
                                                                       20271, 20331, 33451,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 51223, 0, 6,
                                                                       48892, 31066, 48955,
                                                                       20331, 20391, 33541,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 51349, 0, 6,
                                                                       49018, 31246, 49081,
                                                                       20511, 20571, 33811,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 51475, 0, 6,
                                                                       49081, 31291, 49144,
                                                                       20571, 20631, 33901,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 51601, 0, 6,
                                                                       49144, 31336, 49207,
                                                                       20631, 20691, 33991,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 51727, 0, 6,
                                                                       49207, 31381, 49270,
                                                                       20691, 20751, 34081,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 51853, 0, 3, 6,
                                                                       49333, 31741, 49522,
                                                                       50845, 33271, 50971,
                                                                       20871, 21051, 34711,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 52231, 0, 3, 6,
                                                                       49522, 31876, 49711,
                                                                       50971, 33361, 51097,
                                                                       21051, 21231, 34981,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 52609, 0, 3, 6,
                                                                       49711, 32011, 49900,
                                                                       51097, 33451, 51223,
                                                                       21231, 21411, 35251,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 52987, 0, 3, 6,
                                                                       50089, 32551, 50278,
                                                                       51349, 33811, 51475,
                                                                       21771, 21951, 36061,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 53365, 0, 3, 6,
                                                                       50278, 32686, 50467,
                                                                       51475, 33901, 51601,
                                                                       21951, 22131, 36331,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 53743, 0, 3, 6,
                                                                       50467, 32821, 50656,
                                                                       51601, 33991, 51727,
                                                                       22131, 22311, 36601,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54121, 0, 6,
                                                                       50845, 33271, 50971,
                                                                       22671, 22771, 37171,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54331, 0, 6,
                                                                       50971, 33361, 51097,
                                                                       22771, 22871, 37321,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54541, 0, 6,
                                                                       51097, 33451, 51223,
                                                                       22871, 22971, 37471,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54751, 0, 6,
                                                                       51349, 33811, 51475,
                                                                       23171, 23271, 37921,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 54961, 0, 6,
                                                                       51475, 33901, 51601,
                                                                       23271, 23371, 38071,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 55171, 0, 6,
                                                                       51601, 33991, 51727,
                                                                       23371, 23471, 38221,
                                                                       ncols, gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 55381, 0, 3, 6,
                                                                       49333, 49522, 51853,
                                                                       34711, 52231, 54121,
                                                                       37171, 54331, 23671,
                                                                       23971, 39271, ncols,
                                                                       gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 56011, 0, 3, 6,
                                                                       49522, 49711, 52231,
                                                                       34981, 52609, 54331,
                                                                       37321, 54541, 23971,
                                                                       24271, 39721, ncols,
                                                                       gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 56641, 0, 3, 6,
                                                                       50089, 50278, 52987,
                                                                       36061, 53365, 54751,
                                                                       37921, 54961, 24871,
                                                                       25171, 41071, ncols,
                                                                       gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 57271, 0, 3, 6,
                                                                       50278, 50467, 53365,
                                                                       36331, 53743, 54961,
                                                                       38071, 55171, 25171,
                                                                       25471, 41521, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 57901, 0, 6,
                                                                       54121, 37171, 54331,
                                                                       26071, 26221, 42421,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 58216, 0, 6,
                                                                       54331, 37321, 54541,
                                                                       26221, 26371, 42646,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 58531, 0, 6,
                                                                       54751, 37921, 54961,
                                                                       26671, 26821, 43321,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 58846, 0, 6,
                                                                       54961, 38071, 55171,
                                                                       26821, 26971, 43546,
                                                                       ncols, gamma, p, q);

                    compute_prim_gph_three_center_electron_repulsion_0(buffer, 59161, 0, 3, 6,
                                                                       51853, 52231, 55381,
                                                                       39271, 56011, 57901,
                                                                       42421, 58216, 27271,
                                                                       27721, 45121, ncols,
                                                                       gamma, p, q);

                    compute_prim_gph_three_center_electron_repulsion_0(buffer, 60106, 0, 3, 6,
                                                                       52987, 53365, 56641,
                                                                       41071, 57271, 58531,
                                                                       43321, 58846, 28621,
                                                                       29071, 47146, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61051, 6, 29971,
                                                                       29986, 47821, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61079, 6, 29986,
                                                                       30001, 47842, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61107, 6, 30001,
                                                                       30016, 47863, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61135, 6, 30016,
                                                                       30031, 47884, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61163, 6, 30031,
                                                                       30046, 47905, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61191, 6, 30046,
                                                                       30061, 47926, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61219, 6, 30091,
                                                                       30106, 47947, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61247, 6, 30106,
                                                                       30121, 47968, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61275, 6, 30121,
                                                                       30136, 47989, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61303, 6, 30136,
                                                                       30151, 48010, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61331, 6, 30151,
                                                                       30166, 48031, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 61359, 6, 30166,
                                                                       30181, 48052, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61387, 3, 6,
                                                                       61051, 47821, 61079,
                                                                       30211, 30256, 48073,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61471, 3, 6,
                                                                       61079, 47842, 61107,
                                                                       30256, 30301, 48136,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61555, 3, 6,
                                                                       61107, 47863, 61135,
                                                                       30301, 30346, 48199,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61639, 3, 6,
                                                                       61135, 47884, 61163,
                                                                       30346, 30391, 48262,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61723, 3, 6,
                                                                       61163, 47905, 61191,
                                                                       30391, 30436, 48325,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61807, 3, 6,
                                                                       61219, 47947, 61247,
                                                                       30526, 30571, 48388,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61891, 3, 6,
                                                                       61247, 47968, 61275,
                                                                       30571, 30616, 48451,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 61975, 3, 6,
                                                                       61275, 47989, 61303,
                                                                       30616, 30661, 48514,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 62059, 3, 6,
                                                                       61303, 48010, 61331,
                                                                       30661, 30706, 48577,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 62143, 3, 6,
                                                                       61331, 48031, 61359,
                                                                       30706, 30751, 48640,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62227, 0, 6,
                                                                       61051, 47821, 61079,
                                                                       30841, 30886, 48703,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62311, 0, 6,
                                                                       61079, 47842, 61107,
                                                                       30886, 30931, 48766,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62395, 0, 6,
                                                                       61107, 47863, 61135,
                                                                       30931, 30976, 48829,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62479, 0, 6,
                                                                       61135, 47884, 61163,
                                                                       30976, 31021, 48892,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62563, 0, 6,
                                                                       61163, 47905, 61191,
                                                                       31021, 31066, 48955,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62647, 0, 6,
                                                                       61219, 47947, 61247,
                                                                       31156, 31201, 49018,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62731, 0, 6,
                                                                       61247, 47968, 61275,
                                                                       31201, 31246, 49081,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62815, 0, 6,
                                                                       61275, 47989, 61303,
                                                                       31246, 31291, 49144,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62899, 0, 6,
                                                                       61303, 48010, 61331,
                                                                       31291, 31336, 49207,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 62983, 0, 6,
                                                                       61331, 48031, 61359,
                                                                       31336, 31381, 49270,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 63067, 0, 3, 6,
                                                                       61387, 48073, 61471,
                                                                       62227, 48703, 62311,
                                                                       31471, 31606, 49333,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 63319, 0, 3, 6,
                                                                       61471, 48136, 61555,
                                                                       62311, 48766, 62395,
                                                                       31606, 31741, 49522,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 63571, 0, 3, 6,
                                                                       61555, 48199, 61639,
                                                                       62395, 48829, 62479,
                                                                       31741, 31876, 49711,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 63823, 0, 3, 6,
                                                                       61639, 48262, 61723,
                                                                       62479, 48892, 62563,
                                                                       31876, 32011, 49900,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 64075, 0, 3, 6,
                                                                       61807, 48388, 61891,
                                                                       62647, 49018, 62731,
                                                                       32281, 32416, 50089,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 64327, 0, 3, 6,
                                                                       61891, 48451, 61975,
                                                                       62731, 49081, 62815,
                                                                       32416, 32551, 50278,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 64579, 0, 3, 6,
                                                                       61975, 48514, 62059,
                                                                       62815, 49144, 62899,
                                                                       32551, 32686, 50467,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 64831, 0, 3, 6,
                                                                       62059, 48577, 62143,
                                                                       62899, 49207, 62983,
                                                                       32686, 32821, 50656,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 65083, 0, 6,
                                                                       62227, 48703, 62311,
                                                                       33091, 33181, 50845,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 65251, 0, 6,
                                                                       62311, 48766, 62395,
                                                                       33181, 33271, 50971,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 65419, 0, 6,
                                                                       62395, 48829, 62479,
                                                                       33271, 33361, 51097,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 65587, 0, 6,
                                                                       62479, 48892, 62563,
                                                                       33361, 33451, 51223,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 65755, 0, 6,
                                                                       62647, 49018, 62731,
                                                                       33631, 33721, 51349,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 65923, 0, 6,
                                                                       62731, 49081, 62815,
                                                                       33721, 33811, 51475,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 66091, 0, 6,
                                                                       62815, 49144, 62899,
                                                                       33811, 33901, 51601,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 66259, 0, 6,
                                                                       62899, 49207, 62983,
                                                                       33901, 33991, 51727,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 66427, 0, 3, 6,
                                                                       63067, 49333, 63319,
                                                                       65083, 50845, 65251,
                                                                       34171, 34441, 51853,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 66931, 0, 3, 6,
                                                                       63319, 49522, 63571,
                                                                       65251, 50971, 65419,
                                                                       34441, 34711, 52231,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 67435, 0, 3, 6,
                                                                       63571, 49711, 63823,
                                                                       65419, 51097, 65587,
                                                                       34711, 34981, 52609,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 67939, 0, 3, 6,
                                                                       64075, 50089, 64327,
                                                                       65755, 51349, 65923,
                                                                       35521, 35791, 52987,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 68443, 0, 3, 6,
                                                                       64327, 50278, 64579,
                                                                       65923, 51475, 66091,
                                                                       35791, 36061, 53365,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 68947, 0, 3, 6,
                                                                       64579, 50467, 64831,
                                                                       66091, 51601, 66259,
                                                                       36061, 36331, 53743,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 69451, 0, 6,
                                                                       65083, 50845, 65251,
                                                                       36871, 37021, 54121,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 69731, 0, 6,
                                                                       65251, 50971, 65419,
                                                                       37021, 37171, 54331,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 70011, 0, 6,
                                                                       65419, 51097, 65587,
                                                                       37171, 37321, 54541,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 70291, 0, 6,
                                                                       65755, 51349, 65923,
                                                                       37621, 37771, 54751,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 70571, 0, 6,
                                                                       65923, 51475, 66091,
                                                                       37771, 37921, 54961,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 70851, 0, 6,
                                                                       66091, 51601, 66259,
                                                                       37921, 38071, 55171,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpi_three_center_electron_repulsion_0(buffer, 71131, 0, 3, 6,
                                                                       63067, 63319, 66427,
                                                                       51853, 66931, 69451,
                                                                       54121, 69731, 38371,
                                                                       38821, 55381, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpi_three_center_electron_repulsion_0(buffer, 71971, 0, 3, 6,
                                                                       63319, 63571, 66931,
                                                                       52231, 67435, 69731,
                                                                       54331, 70011, 38821,
                                                                       39271, 56011, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpi_three_center_electron_repulsion_0(buffer, 72811, 0, 3, 6,
                                                                       64075, 64327, 67939,
                                                                       52987, 68443, 70291,
                                                                       54751, 70571, 40171,
                                                                       40621, 56641, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpi_three_center_electron_repulsion_0(buffer, 73651, 0, 3, 6,
                                                                       64327, 64579, 68443,
                                                                       53365, 68947, 70571,
                                                                       54961, 70851, 40621,
                                                                       41071, 57271, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 74491, 0, 6,
                                                                       69451, 54121, 69731,
                                                                       41971, 42196, 57901,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 74911, 0, 6,
                                                                       69731, 54331, 70011,
                                                                       42196, 42421, 58216,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 75331, 0, 6,
                                                                       70291, 54751, 70571,
                                                                       42871, 43096, 58531,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 75751, 0, 6,
                                                                       70571, 54961, 70851,
                                                                       43096, 43321, 58846,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpi_three_center_electron_repulsion_0(buffer, 76171, 0, 3, 6,
                                                                       66427, 66931, 71131,
                                                                       55381, 71971, 74491,
                                                                       57901, 74911, 43771,
                                                                       44446, 59161, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpi_three_center_electron_repulsion_0(buffer, 77431, 0, 3, 6,
                                                                       67939, 68443, 72811,
                                                                       56641, 73651, 75331,
                                                                       58531, 75751, 45796,
                                                                       46471, 60106, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 78691, 77431, 15, 28, ncols, beta);

                    simdgeo::geom_s_y(buffer, 79111, 77431, 15, 28, ncols, beta);

                    simdgeo::geom_s_z(buffer, 79531, 77431, 15, 28, ncols, beta);

                    simdgeo::geom_s_x(buffer, 79951, 76171, 15, 28, ncols, beta);

                    simdgeo::geom_s_y(buffer, 80371, 76171, 15, 28, ncols, beta);

                    simdgeo::geom_s_z(buffer, 80791, 76171, 15, 28, ncols, beta);

                    simdfunc::contract_primitives(buffer, 81211, 78691, 2520, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 83731, 81211, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 83731, 13, nmax);

        simdtrf::transform_i_inner(buffer, 83731, 81631, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 117 * nvalues + n * npairs, nvalues, buffer, 83731,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 83731, 82051, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 234 * nvalues + n * npairs, nvalues, buffer, 83731,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 83731, 82471, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 351 * nvalues + n * npairs, nvalues, buffer, 83731,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 83731, 82891, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 468 * nvalues + n * npairs, nvalues, buffer, 83731,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 83731, 83311, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 585 * nvalues + n * npairs, nvalues, buffer, 83731,
                                   13, nmax);
    }

    for (size_t m = 0; m < 702; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
