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


#include "SimdThreeCenterElectronRepulsionGeom010RecGSI.hpp"

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
compute_geom_010_gsi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_geom_010_gsi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 42065, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 351 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 42065, 40610, 1260, dimensions);

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

                    simdfunc::compute_full_t3c_boys_function(buffer, coordinates, 9, 6, 11,
                                                             ncols, fj, i * nprim_b + j, fq);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 22, 3, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 25, 3, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 28, 3, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 31, 3, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 34, 3, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 37, 3, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 40, 3, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 43, 3, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 46, 3, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 49, 3, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_sps_three_center_electron_repulsion_0(buffer, 52, 3, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 55, 0, 6, 10, 11,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 58, 0, 6, 11, 12,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 61, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 64, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 67, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 70, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 73, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 76, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 79, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 82, 0, 6, 19, 20,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 85, 0, 6, 20, 21,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 88, 0, 6, 10, 11,
                                                                       22, 25, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 97, 0, 6, 11, 12,
                                                                       25, 28, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 106, 0, 6, 12, 13,
                                                                       28, 31, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 115, 0, 6, 13, 14,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 124, 0, 6, 14, 15,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 133, 0, 6, 15, 16,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 142, 0, 6, 16, 17,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 151, 0, 6, 17, 18,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 160, 0, 6, 18, 19,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 169, 0, 6, 19, 20,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 178, 0, 6, 10, 11,
                                                                       55, 58, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 184, 0, 6, 11, 12,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 190, 0, 6, 12, 13,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 196, 0, 6, 13, 14,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 202, 0, 6, 14, 15,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 208, 0, 6, 15, 16,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 214, 0, 6, 16, 17,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 220, 0, 6, 17, 18,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 226, 0, 6, 18, 19,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_dss_three_center_electron_repulsion_0(buffer, 232, 0, 6, 19, 20,
                                                                       82, 85, ncols, gamma, p,
                                                                       q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 238, 0, 3, 6, 55,
                                                                       58, 88, 97, 178, 184,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 256, 0, 3, 6, 58,
                                                                       61, 97, 106, 184, 190,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 274, 0, 3, 6, 61,
                                                                       64, 106, 115, 190, 196,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 292, 0, 3, 6, 64,
                                                                       67, 115, 124, 196, 202,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 310, 0, 3, 6, 67,
                                                                       70, 124, 133, 202, 208,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 328, 0, 3, 6, 70,
                                                                       73, 133, 142, 208, 214,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 346, 0, 3, 6, 73,
                                                                       76, 142, 151, 214, 220,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 364, 0, 3, 6, 76,
                                                                       79, 151, 160, 220, 226,
                                                                       ncols, gamma, p, q);

                    compute_prim_dps_three_center_electron_repulsion_0(buffer, 382, 0, 3, 6, 79,
                                                                       82, 160, 169, 226, 232,
                                                                       ncols, gamma, p, q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 400, 0, 6, 55, 58,
                                                                       178, 184, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 410, 0, 6, 58, 61,
                                                                       184, 190, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 420, 0, 6, 61, 64,
                                                                       190, 196, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 430, 0, 6, 64, 67,
                                                                       196, 202, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 440, 0, 6, 67, 70,
                                                                       202, 208, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 450, 0, 6, 70, 73,
                                                                       208, 214, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 460, 0, 6, 73, 76,
                                                                       214, 220, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 470, 0, 6, 76, 79,
                                                                       220, 226, ncols, gamma, p,
                                                                       q);

                    compute_prim_fss_three_center_electron_repulsion_0(buffer, 480, 0, 6, 79, 82,
                                                                       226, 232, ncols, gamma, p,
                                                                       q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 490, 0, 3, 6, 178,
                                                                       184, 238, 256, 400, 410,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 520, 0, 3, 6, 184,
                                                                       190, 256, 274, 410, 420,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 550, 0, 3, 6, 190,
                                                                       196, 274, 292, 420, 430,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 580, 0, 3, 6, 196,
                                                                       202, 292, 310, 430, 440,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 610, 0, 3, 6, 202,
                                                                       208, 310, 328, 440, 450,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 640, 0, 3, 6, 208,
                                                                       214, 328, 346, 450, 460,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 670, 0, 3, 6, 214,
                                                                       220, 346, 364, 460, 470,
                                                                       ncols, gamma, p, q);

                    compute_prim_fps_three_center_electron_repulsion_0(buffer, 700, 0, 3, 6, 220,
                                                                       226, 364, 382, 470, 480,
                                                                       ncols, gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 730, 0, 6, 178,
                                                                       184, 400, 410, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 745, 0, 6, 184,
                                                                       190, 410, 420, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 760, 0, 6, 190,
                                                                       196, 420, 430, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 775, 0, 6, 196,
                                                                       202, 430, 440, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 790, 0, 6, 202,
                                                                       208, 440, 450, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 805, 0, 6, 208,
                                                                       214, 450, 460, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 820, 0, 6, 214,
                                                                       220, 460, 470, ncols,
                                                                       gamma, p, q);

                    compute_prim_gss_three_center_electron_repulsion_0(buffer, 835, 0, 6, 220,
                                                                       226, 470, 480, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 850, 0, 3, 6, 238,
                                                                       256, 400, 410, 490, 520,
                                                                       730, 745, ncols, gamma, p,
                                                                       q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 895, 0, 3, 6, 256,
                                                                       274, 410, 420, 520, 550,
                                                                       745, 760, ncols, gamma, p,
                                                                       q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 940, 0, 3, 6, 274,
                                                                       292, 420, 430, 550, 580,
                                                                       760, 775, ncols, gamma, p,
                                                                       q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 985, 0, 3, 6, 292,
                                                                       310, 430, 440, 580, 610,
                                                                       775, 790, ncols, gamma, p,
                                                                       q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1030, 0, 3, 6,
                                                                       310, 328, 440, 450, 610,
                                                                       640, 790, 805, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1075, 0, 3, 6,
                                                                       328, 346, 450, 460, 640,
                                                                       670, 805, 820, ncols,
                                                                       gamma, p, q);

                    compute_prim_gps_three_center_electron_repulsion_0(buffer, 1120, 0, 3, 6,
                                                                       346, 364, 460, 470, 670,
                                                                       700, 820, 835, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1165, 6, 12,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1168, 6, 13,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1171, 6, 14,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1174, 6, 15,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1177, 6, 16,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1180, 6, 17,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1183, 6, 18,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1186, 6, 19,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1189, 6, 20,
                                                                       ncols, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 1192, 6, 21,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1195, 6, 12, 28,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1204, 6, 13, 31,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1213, 6, 14, 34,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1222, 6, 15, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1231, 6, 16, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1240, 6, 17, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1249, 6, 18, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1258, 6, 19, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 1267, 6, 20, 52,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1276, 6, 12, 61,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1285, 6, 13, 64,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1294, 6, 14, 67,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1303, 6, 15, 70,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1312, 6, 16, 73,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1321, 6, 17, 76,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1330, 6, 18, 79,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1339, 6, 19, 82,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1348, 6, 20, 85,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1357, 6, 28, 61,
                                                                       106, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1384, 6, 31, 64,
                                                                       115, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1411, 6, 34, 67,
                                                                       124, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1438, 6, 37, 70,
                                                                       133, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1465, 6, 40, 73,
                                                                       142, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1492, 6, 43, 76,
                                                                       151, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1519, 6, 46, 79,
                                                                       160, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1546, 6, 49, 82,
                                                                       169, ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1573, 6, 61, 190,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1591, 6, 64, 196,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1609, 6, 67, 202,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1627, 6, 70, 208,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1645, 6, 73, 214,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1663, 6, 76, 220,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1681, 6, 79, 226,
                                                                       ncols, p, q);

                    compute_prim_dsp_three_center_electron_repulsion_0(buffer, 1699, 6, 82, 232,
                                                                       ncols, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1717, 0, 6, 1357,
                                                                       106, 1384, 190, 274,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1771, 0, 6, 1384,
                                                                       115, 1411, 196, 292,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1825, 0, 6, 1411,
                                                                       124, 1438, 202, 310,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1879, 0, 6, 1438,
                                                                       133, 1465, 208, 328,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1933, 0, 6, 1465,
                                                                       142, 1492, 214, 346,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 1987, 0, 6, 1492,
                                                                       151, 1519, 220, 364,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpp_three_center_electron_repulsion_0(buffer, 2041, 0, 6, 1519,
                                                                       160, 1546, 226, 382,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2095, 6, 190, 420,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2125, 6, 196, 430,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2155, 6, 202, 440,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2185, 6, 208, 450,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2215, 6, 214, 460,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2245, 6, 220, 470,
                                                                       ncols, p, q);

                    compute_prim_fsp_three_center_electron_repulsion_0(buffer, 2275, 6, 226, 480,
                                                                       ncols, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2305, 0, 6, 1717,
                                                                       274, 1771, 420, 550,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2395, 0, 6, 1771,
                                                                       292, 1825, 430, 580,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2485, 0, 6, 1825,
                                                                       310, 1879, 440, 610,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2575, 0, 6, 1879,
                                                                       328, 1933, 450, 640,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2665, 0, 6, 1933,
                                                                       346, 1987, 460, 670,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpp_three_center_electron_repulsion_0(buffer, 2755, 0, 6, 1987,
                                                                       364, 2041, 470, 700,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2845, 6, 420, 760,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2890, 6, 430, 775,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2935, 6, 440, 790,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 2980, 6, 450, 805,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3025, 6, 460, 820,
                                                                       ncols, p, q);

                    compute_prim_gsp_three_center_electron_repulsion_0(buffer, 3070, 6, 470, 835,
                                                                       ncols, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 3115, 0, 6, 2305,
                                                                       550, 2395, 760, 940,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 3250, 0, 6, 2395,
                                                                       580, 2485, 775, 985,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 3385, 0, 6, 2485,
                                                                       610, 2575, 790, 1030,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 3520, 0, 6, 2575,
                                                                       640, 2665, 805, 1075,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpp_three_center_electron_repulsion_0(buffer, 3655, 0, 6, 2665,
                                                                       670, 2755, 820, 1120,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3790, 6, 10, 11,
                                                                       1165, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3796, 6, 11, 12,
                                                                       1168, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3802, 6, 12, 13,
                                                                       1171, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3808, 6, 13, 14,
                                                                       1174, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3814, 6, 14, 15,
                                                                       1177, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3820, 6, 15, 16,
                                                                       1180, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3826, 6, 16, 17,
                                                                       1183, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3832, 6, 17, 18,
                                                                       1186, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3838, 6, 18, 19,
                                                                       1189, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 3844, 6, 19, 20,
                                                                       1192, ncols, gamma, p,
                                                                       q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3850, 3, 6, 3790,
                                                                       1165, 3796, 1195, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3868, 3, 6, 3796,
                                                                       1168, 3802, 1204, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3886, 3, 6, 3802,
                                                                       1171, 3808, 1213, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3904, 3, 6, 3808,
                                                                       1174, 3814, 1222, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3922, 3, 6, 3814,
                                                                       1177, 3820, 1231, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3940, 3, 6, 3820,
                                                                       1180, 3826, 1240, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3958, 3, 6, 3826,
                                                                       1183, 3832, 1249, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3976, 3, 6, 3832,
                                                                       1186, 3838, 1258, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 3994, 3, 6, 3838,
                                                                       1189, 3844, 1267, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4012, 0, 6, 3790,
                                                                       1165, 3796, 1276, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4030, 0, 6, 3796,
                                                                       1168, 3802, 1285, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4048, 0, 6, 3802,
                                                                       1171, 3808, 1294, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4066, 0, 6, 3808,
                                                                       1174, 3814, 1303, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4084, 0, 6, 3814,
                                                                       1177, 3820, 1312, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4102, 0, 6, 3820,
                                                                       1180, 3826, 1321, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4120, 0, 6, 3826,
                                                                       1183, 3832, 1330, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4138, 0, 6, 3832,
                                                                       1186, 3838, 1339, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 4156, 0, 6, 3838,
                                                                       1189, 3844, 1348, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4174, 0, 3, 6,
                                                                       3850, 1195, 3868, 4012,
                                                                       1276, 4030, 88, 97, 1357,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4228, 0, 3, 6,
                                                                       3868, 1204, 3886, 4030,
                                                                       1285, 4048, 97, 106, 1384,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4282, 0, 3, 6,
                                                                       3886, 1213, 3904, 4048,
                                                                       1294, 4066, 106, 115,
                                                                       1411, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4336, 0, 3, 6,
                                                                       3904, 1222, 3922, 4066,
                                                                       1303, 4084, 115, 124,
                                                                       1438, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4390, 0, 3, 6,
                                                                       3922, 1231, 3940, 4084,
                                                                       1312, 4102, 124, 133,
                                                                       1465, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4444, 0, 3, 6,
                                                                       3940, 1240, 3958, 4102,
                                                                       1321, 4120, 133, 142,
                                                                       1492, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4498, 0, 3, 6,
                                                                       3958, 1249, 3976, 4120,
                                                                       1330, 4138, 142, 151,
                                                                       1519, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 4552, 0, 3, 6,
                                                                       3976, 1258, 3994, 4138,
                                                                       1339, 4156, 151, 160,
                                                                       1546, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4606, 0, 6, 4012,
                                                                       1276, 4030, 178, 184,
                                                                       1573, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4642, 0, 6, 4030,
                                                                       1285, 4048, 184, 190,
                                                                       1591, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4678, 0, 6, 4048,
                                                                       1294, 4066, 190, 196,
                                                                       1609, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4714, 0, 6, 4066,
                                                                       1303, 4084, 196, 202,
                                                                       1627, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4750, 0, 6, 4084,
                                                                       1312, 4102, 202, 208,
                                                                       1645, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4786, 0, 6, 4102,
                                                                       1321, 4120, 208, 214,
                                                                       1663, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4822, 0, 6, 4120,
                                                                       1330, 4138, 214, 220,
                                                                       1681, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsd_three_center_electron_repulsion_0(buffer, 4858, 0, 6, 4138,
                                                                       1339, 4156, 220, 226,
                                                                       1699, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 4894, 0, 3, 6,
                                                                       4174, 1357, 4228, 4606,
                                                                       1573, 4642, 238, 256,
                                                                       1717, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5002, 0, 3, 6,
                                                                       4228, 1384, 4282, 4642,
                                                                       1591, 4678, 256, 274,
                                                                       1771, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5110, 0, 3, 6,
                                                                       4282, 1411, 4336, 4678,
                                                                       1609, 4714, 274, 292,
                                                                       1825, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5218, 0, 3, 6,
                                                                       4336, 1438, 4390, 4714,
                                                                       1627, 4750, 292, 310,
                                                                       1879, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5326, 0, 3, 6,
                                                                       4390, 1465, 4444, 4750,
                                                                       1645, 4786, 310, 328,
                                                                       1933, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5434, 0, 3, 6,
                                                                       4444, 1492, 4498, 4786,
                                                                       1663, 4822, 328, 346,
                                                                       1987, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpd_three_center_electron_repulsion_0(buffer, 5542, 0, 3, 6,
                                                                       4498, 1519, 4552, 4822,
                                                                       1681, 4858, 346, 364,
                                                                       2041, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5650, 0, 6, 4606,
                                                                       1573, 4642, 400, 410,
                                                                       2095, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5710, 0, 6, 4642,
                                                                       1591, 4678, 410, 420,
                                                                       2125, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5770, 0, 6, 4678,
                                                                       1609, 4714, 420, 430,
                                                                       2155, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5830, 0, 6, 4714,
                                                                       1627, 4750, 430, 440,
                                                                       2185, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5890, 0, 6, 4750,
                                                                       1645, 4786, 440, 450,
                                                                       2215, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 5950, 0, 6, 4786,
                                                                       1663, 4822, 450, 460,
                                                                       2245, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsd_three_center_electron_repulsion_0(buffer, 6010, 0, 6, 4822,
                                                                       1681, 4858, 460, 470,
                                                                       2275, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6070, 0, 3, 6,
                                                                       4174, 4228, 4894, 1717,
                                                                       5002, 5650, 2095, 5710,
                                                                       490, 520, 2305, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6250, 0, 3, 6,
                                                                       4228, 4282, 5002, 1771,
                                                                       5110, 5710, 2125, 5770,
                                                                       520, 550, 2395, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6430, 0, 3, 6,
                                                                       4282, 4336, 5110, 1825,
                                                                       5218, 5770, 2155, 5830,
                                                                       550, 580, 2485, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6610, 0, 3, 6,
                                                                       4336, 4390, 5218, 1879,
                                                                       5326, 5830, 2185, 5890,
                                                                       580, 610, 2575, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6790, 0, 3, 6,
                                                                       4390, 4444, 5326, 1933,
                                                                       5434, 5890, 2215, 5950,
                                                                       610, 640, 2665, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpd_three_center_electron_repulsion_0(buffer, 6970, 0, 3, 6,
                                                                       4444, 4498, 5434, 1987,
                                                                       5542, 5950, 2245, 6010,
                                                                       640, 670, 2755, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7150, 0, 6, 5650,
                                                                       2095, 5710, 730, 745,
                                                                       2845, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7240, 0, 6, 5710,
                                                                       2125, 5770, 745, 760,
                                                                       2890, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7330, 0, 6, 5770,
                                                                       2155, 5830, 760, 775,
                                                                       2935, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7420, 0, 6, 5830,
                                                                       2185, 5890, 775, 790,
                                                                       2980, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7510, 0, 6, 5890,
                                                                       2215, 5950, 790, 805,
                                                                       3025, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsd_three_center_electron_repulsion_0(buffer, 7600, 0, 6, 5950,
                                                                       2245, 6010, 805, 820,
                                                                       3070, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 7690, 0, 3, 6,
                                                                       4894, 5002, 6070, 2305,
                                                                       6250, 7150, 2845, 7240,
                                                                       850, 895, 3115, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 7960, 0, 3, 6,
                                                                       5002, 5110, 6250, 2395,
                                                                       6430, 7240, 2890, 7330,
                                                                       895, 940, 3250, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 8230, 0, 3, 6,
                                                                       5110, 5218, 6430, 2485,
                                                                       6610, 7330, 2935, 7420,
                                                                       940, 985, 3385, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 8500, 0, 3, 6,
                                                                       5218, 5326, 6610, 2575,
                                                                       6790, 7420, 2980, 7510,
                                                                       985, 1030, 3520, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpd_three_center_electron_repulsion_0(buffer, 8770, 0, 3, 6,
                                                                       5326, 5434, 6790, 2665,
                                                                       6970, 7510, 3025, 7600,
                                                                       1030, 1075, 3655, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9040, 6, 1165,
                                                                       1168, 3802, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9050, 6, 1168,
                                                                       1171, 3808, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9060, 6, 1171,
                                                                       1174, 3814, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9070, 6, 1174,
                                                                       1177, 3820, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9080, 6, 1177,
                                                                       1180, 3826, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9090, 6, 1180,
                                                                       1183, 3832, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9100, 6, 1183,
                                                                       1186, 3838, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 9110, 6, 1186,
                                                                       1189, 3844, ncols, gamma,
                                                                       p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9120, 3, 6, 9040,
                                                                       3802, 9050, 3886, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9150, 3, 6, 9050,
                                                                       3808, 9060, 3904, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9180, 3, 6, 9060,
                                                                       3814, 9070, 3922, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9210, 3, 6, 9070,
                                                                       3820, 9080, 3940, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9240, 3, 6, 9080,
                                                                       3826, 9090, 3958, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9270, 3, 6, 9090,
                                                                       3832, 9100, 3976, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 9300, 3, 6, 9100,
                                                                       3838, 9110, 3994, ncols,
                                                                       gamma, p, q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9330, 0, 6, 9040,
                                                                       3802, 9050, 1276, 1285,
                                                                       4048, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9360, 0, 6, 9050,
                                                                       3808, 9060, 1285, 1294,
                                                                       4066, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9390, 0, 6, 9060,
                                                                       3814, 9070, 1294, 1303,
                                                                       4084, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9420, 0, 6, 9070,
                                                                       3820, 9080, 1303, 1312,
                                                                       4102, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9450, 0, 6, 9080,
                                                                       3826, 9090, 1312, 1321,
                                                                       4120, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9480, 0, 6, 9090,
                                                                       3832, 9100, 1321, 1330,
                                                                       4138, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 9510, 0, 6, 9100,
                                                                       3838, 9110, 1330, 1339,
                                                                       4156, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9540, 0, 3, 6,
                                                                       9120, 3886, 9150, 9330,
                                                                       4048, 9360, 1357, 1384,
                                                                       4282, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9630, 0, 3, 6,
                                                                       9150, 3904, 9180, 9360,
                                                                       4066, 9390, 1384, 1411,
                                                                       4336, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9720, 0, 3, 6,
                                                                       9180, 3922, 9210, 9390,
                                                                       4084, 9420, 1411, 1438,
                                                                       4390, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9810, 0, 3, 6,
                                                                       9210, 3940, 9240, 9420,
                                                                       4102, 9450, 1438, 1465,
                                                                       4444, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9900, 0, 3, 6,
                                                                       9240, 3958, 9270, 9450,
                                                                       4120, 9480, 1465, 1492,
                                                                       4498, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 9990, 0, 3, 6,
                                                                       9270, 3976, 9300, 9480,
                                                                       4138, 9510, 1492, 1519,
                                                                       4552, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10080, 0, 6, 9330,
                                                                       4048, 9360, 1573, 1591,
                                                                       4678, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10140, 0, 6, 9360,
                                                                       4066, 9390, 1591, 1609,
                                                                       4714, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10200, 0, 6, 9390,
                                                                       4084, 9420, 1609, 1627,
                                                                       4750, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10260, 0, 6, 9420,
                                                                       4102, 9450, 1627, 1645,
                                                                       4786, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10320, 0, 6, 9450,
                                                                       4120, 9480, 1645, 1663,
                                                                       4822, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsf_three_center_electron_repulsion_0(buffer, 10380, 0, 6, 9480,
                                                                       4138, 9510, 1663, 1681,
                                                                       4858, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10440, 0, 3, 6,
                                                                       9540, 4282, 9630, 10080,
                                                                       4678, 10140, 1717, 1771,
                                                                       5110, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10620, 0, 3, 6,
                                                                       9630, 4336, 9720, 10140,
                                                                       4714, 10200, 1771, 1825,
                                                                       5218, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10800, 0, 3, 6,
                                                                       9720, 4390, 9810, 10200,
                                                                       4750, 10260, 1825, 1879,
                                                                       5326, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 10980, 0, 3, 6,
                                                                       9810, 4444, 9900, 10260,
                                                                       4786, 10320, 1879, 1933,
                                                                       5434, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpf_three_center_electron_repulsion_0(buffer, 11160, 0, 3, 6,
                                                                       9900, 4498, 9990, 10320,
                                                                       4822, 10380, 1933, 1987,
                                                                       5542, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11340, 0, 6,
                                                                       10080, 4678, 10140, 2095,
                                                                       2125, 5770, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11440, 0, 6,
                                                                       10140, 4714, 10200, 2125,
                                                                       2155, 5830, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11540, 0, 6,
                                                                       10200, 4750, 10260, 2155,
                                                                       2185, 5890, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11640, 0, 6,
                                                                       10260, 4786, 10320, 2185,
                                                                       2215, 5950, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsf_three_center_electron_repulsion_0(buffer, 11740, 0, 6,
                                                                       10320, 4822, 10380, 2215,
                                                                       2245, 6010, ncols, gamma,
                                                                       p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 11840, 0, 3, 6,
                                                                       9540, 9630, 10440, 5110,
                                                                       10620, 11340, 5770, 11440,
                                                                       2305, 2395, 6430, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 12140, 0, 3, 6,
                                                                       9630, 9720, 10620, 5218,
                                                                       10800, 11440, 5830, 11540,
                                                                       2395, 2485, 6610, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 12440, 0, 3, 6,
                                                                       9720, 9810, 10800, 5326,
                                                                       10980, 11540, 5890, 11640,
                                                                       2485, 2575, 6790, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpf_three_center_electron_repulsion_0(buffer, 12740, 0, 3, 6,
                                                                       9810, 9900, 10980, 5434,
                                                                       11160, 11640, 5950, 11740,
                                                                       2575, 2665, 6970, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13040, 0, 6,
                                                                       11340, 5770, 11440, 2845,
                                                                       2890, 7330, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13190, 0, 6,
                                                                       11440, 5830, 11540, 2890,
                                                                       2935, 7420, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13340, 0, 6,
                                                                       11540, 5890, 11640, 2935,
                                                                       2980, 7510, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsf_three_center_electron_repulsion_0(buffer, 13490, 0, 6,
                                                                       11640, 5950, 11740, 2980,
                                                                       3025, 7600, ncols, gamma,
                                                                       p, q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 13640, 0, 3, 6,
                                                                       10440, 10620, 11840, 6430,
                                                                       12140, 13040, 7330, 13190,
                                                                       3115, 3250, 8230, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 14090, 0, 3, 6,
                                                                       10620, 10800, 12140, 6610,
                                                                       12440, 13190, 7420, 13340,
                                                                       3250, 3385, 8500, ncols,
                                                                       gamma, p, q);

                    compute_prim_gpf_three_center_electron_repulsion_0(buffer, 14540, 0, 3, 6,
                                                                       10800, 10980, 12440, 6790,
                                                                       12740, 13340, 7510, 13490,
                                                                       3385, 3520, 8770, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 14990, 6, 3790,
                                                                       3796, 9040, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15005, 6, 3796,
                                                                       3802, 9050, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15020, 6, 3802,
                                                                       3808, 9060, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15035, 6, 3808,
                                                                       3814, 9070, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15050, 6, 3814,
                                                                       3820, 9080, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15065, 6, 3820,
                                                                       3826, 9090, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15080, 6, 3826,
                                                                       3832, 9100, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 15095, 6, 3832,
                                                                       3838, 9110, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15110, 3, 6,
                                                                       14990, 9040, 15005, 3850,
                                                                       3868, 9120, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15155, 3, 6,
                                                                       15005, 9050, 15020, 3868,
                                                                       3886, 9150, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15200, 3, 6,
                                                                       15020, 9060, 15035, 3886,
                                                                       3904, 9180, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15245, 3, 6,
                                                                       15035, 9070, 15050, 3904,
                                                                       3922, 9210, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15290, 3, 6,
                                                                       15050, 9080, 15065, 3922,
                                                                       3940, 9240, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15335, 3, 6,
                                                                       15065, 9090, 15080, 3940,
                                                                       3958, 9270, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 15380, 3, 6,
                                                                       15080, 9100, 15095, 3958,
                                                                       3976, 9300, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15425, 0, 6,
                                                                       14990, 9040, 15005, 4012,
                                                                       4030, 9330, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15470, 0, 6,
                                                                       15005, 9050, 15020, 4030,
                                                                       4048, 9360, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15515, 0, 6,
                                                                       15020, 9060, 15035, 4048,
                                                                       4066, 9390, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15560, 0, 6,
                                                                       15035, 9070, 15050, 4066,
                                                                       4084, 9420, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15605, 0, 6,
                                                                       15050, 9080, 15065, 4084,
                                                                       4102, 9450, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15650, 0, 6,
                                                                       15065, 9090, 15080, 4102,
                                                                       4120, 9480, ncols, gamma,
                                                                       p, q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 15695, 0, 6,
                                                                       15080, 9100, 15095, 4120,
                                                                       4138, 9510, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 15740, 0, 3, 6,
                                                                       15110, 9120, 15155, 15425,
                                                                       9330, 15470, 4174, 4228,
                                                                       9540, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 15875, 0, 3, 6,
                                                                       15155, 9150, 15200, 15470,
                                                                       9360, 15515, 4228, 4282,
                                                                       9630, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 16010, 0, 3, 6,
                                                                       15200, 9180, 15245, 15515,
                                                                       9390, 15560, 4282, 4336,
                                                                       9720, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 16145, 0, 3, 6,
                                                                       15245, 9210, 15290, 15560,
                                                                       9420, 15605, 4336, 4390,
                                                                       9810, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 16280, 0, 3, 6,
                                                                       15290, 9240, 15335, 15605,
                                                                       9450, 15650, 4390, 4444,
                                                                       9900, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 16415, 0, 3, 6,
                                                                       15335, 9270, 15380, 15650,
                                                                       9480, 15695, 4444, 4498,
                                                                       9990, ncols, gamma, p,
                                                                       q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 16550, 0, 6,
                                                                       15425, 9330, 15470, 4606,
                                                                       4642, 10080, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 16640, 0, 6,
                                                                       15470, 9360, 15515, 4642,
                                                                       4678, 10140, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 16730, 0, 6,
                                                                       15515, 9390, 15560, 4678,
                                                                       4714, 10200, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 16820, 0, 6,
                                                                       15560, 9420, 15605, 4714,
                                                                       4750, 10260, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 16910, 0, 6,
                                                                       15605, 9450, 15650, 4750,
                                                                       4786, 10320, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsg_three_center_electron_repulsion_0(buffer, 17000, 0, 6,
                                                                       15650, 9480, 15695, 4786,
                                                                       4822, 10380, ncols, gamma,
                                                                       p, q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 17090, 0, 3, 6,
                                                                       15740, 9540, 15875, 16550,
                                                                       10080, 16640, 4894, 5002,
                                                                       10440, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 17360, 0, 3, 6,
                                                                       15875, 9630, 16010, 16640,
                                                                       10140, 16730, 5002, 5110,
                                                                       10620, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 17630, 0, 3, 6,
                                                                       16010, 9720, 16145, 16730,
                                                                       10200, 16820, 5110, 5218,
                                                                       10800, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 17900, 0, 3, 6,
                                                                       16145, 9810, 16280, 16820,
                                                                       10260, 16910, 5218, 5326,
                                                                       10980, ncols, gamma, p,
                                                                       q);

                    compute_prim_dpg_three_center_electron_repulsion_0(buffer, 18170, 0, 3, 6,
                                                                       16280, 9900, 16415, 16910,
                                                                       10320, 17000, 5326, 5434,
                                                                       11160, ncols, gamma, p,
                                                                       q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18440, 0, 6,
                                                                       16550, 10080, 16640, 5650,
                                                                       5710, 11340, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18590, 0, 6,
                                                                       16640, 10140, 16730, 5710,
                                                                       5770, 11440, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18740, 0, 6,
                                                                       16730, 10200, 16820, 5770,
                                                                       5830, 11540, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 18890, 0, 6,
                                                                       16820, 10260, 16910, 5830,
                                                                       5890, 11640, ncols, gamma,
                                                                       p, q);

                    compute_prim_fsg_three_center_electron_repulsion_0(buffer, 19040, 0, 6,
                                                                       16910, 10320, 17000, 5890,
                                                                       5950, 11740, ncols, gamma,
                                                                       p, q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 19190, 0, 3, 6,
                                                                       15740, 15875, 17090,
                                                                       10440, 17360, 18440,
                                                                       11340, 18590, 6070, 6250,
                                                                       11840, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 19640, 0, 3, 6,
                                                                       15875, 16010, 17360,
                                                                       10620, 17630, 18590,
                                                                       11440, 18740, 6250, 6430,
                                                                       12140, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 20090, 0, 3, 6,
                                                                       16010, 16145, 17630,
                                                                       10800, 17900, 18740,
                                                                       11540, 18890, 6430, 6610,
                                                                       12440, ncols, gamma, p,
                                                                       q);

                    compute_prim_fpg_three_center_electron_repulsion_0(buffer, 20540, 0, 3, 6,
                                                                       16145, 16280, 17900,
                                                                       10980, 18170, 18890,
                                                                       11640, 19040, 6610, 6790,
                                                                       12740, ncols, gamma, p,
                                                                       q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 20990, 0, 6,
                                                                       18440, 11340, 18590, 7150,
                                                                       7240, 13040, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 21215, 0, 6,
                                                                       18590, 11440, 18740, 7240,
                                                                       7330, 13190, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 21440, 0, 6,
                                                                       18740, 11540, 18890, 7330,
                                                                       7420, 13340, ncols, gamma,
                                                                       p, q);

                    compute_prim_gsg_three_center_electron_repulsion_0(buffer, 21665, 0, 6,
                                                                       18890, 11640, 19040, 7420,
                                                                       7510, 13490, ncols, gamma,
                                                                       p, q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 21890, 0, 3, 6,
                                                                       17090, 17360, 19190,
                                                                       11840, 19640, 20990,
                                                                       13040, 21215, 7690, 7960,
                                                                       13640, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 22565, 0, 3, 6,
                                                                       17360, 17630, 19640,
                                                                       12140, 20090, 21215,
                                                                       13190, 21440, 7960, 8230,
                                                                       14090, ncols, gamma, p,
                                                                       q);

                    compute_prim_gpg_three_center_electron_repulsion_0(buffer, 23240, 0, 3, 6,
                                                                       17630, 17900, 20090,
                                                                       12440, 20540, 21440,
                                                                       13340, 21665, 8230, 8500,
                                                                       14540, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 23915, 6, 9040,
                                                                       9050, 15020, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 23936, 6, 9050,
                                                                       9060, 15035, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 23957, 6, 9060,
                                                                       9070, 15050, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 23978, 6, 9070,
                                                                       9080, 15065, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 23999, 6, 9080,
                                                                       9090, 15080, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 24020, 6, 9090,
                                                                       9100, 15095, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 24041, 3, 6,
                                                                       23915, 15020, 23936, 9120,
                                                                       9150, 15200, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 24104, 3, 6,
                                                                       23936, 15035, 23957, 9150,
                                                                       9180, 15245, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 24167, 3, 6,
                                                                       23957, 15050, 23978, 9180,
                                                                       9210, 15290, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 24230, 3, 6,
                                                                       23978, 15065, 23999, 9210,
                                                                       9240, 15335, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 24293, 3, 6,
                                                                       23999, 15080, 24020, 9240,
                                                                       9270, 15380, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 24356, 0, 6,
                                                                       23915, 15020, 23936, 9330,
                                                                       9360, 15515, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 24419, 0, 6,
                                                                       23936, 15035, 23957, 9360,
                                                                       9390, 15560, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 24482, 0, 6,
                                                                       23957, 15050, 23978, 9390,
                                                                       9420, 15605, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 24545, 0, 6,
                                                                       23978, 15065, 23999, 9420,
                                                                       9450, 15650, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 24608, 0, 6,
                                                                       23999, 15080, 24020, 9450,
                                                                       9480, 15695, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 24671, 0, 3, 6,
                                                                       24041, 15200, 24104,
                                                                       24356, 15515, 24419, 9540,
                                                                       9630, 16010, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 24860, 0, 3, 6,
                                                                       24104, 15245, 24167,
                                                                       24419, 15560, 24482, 9630,
                                                                       9720, 16145, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 25049, 0, 3, 6,
                                                                       24167, 15290, 24230,
                                                                       24482, 15605, 24545, 9720,
                                                                       9810, 16280, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 25238, 0, 3, 6,
                                                                       24230, 15335, 24293,
                                                                       24545, 15650, 24608, 9810,
                                                                       9900, 16415, ncols, gamma,
                                                                       p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 25427, 0, 6,
                                                                       24356, 15515, 24419,
                                                                       10080, 10140, 16730,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 25553, 0, 6,
                                                                       24419, 15560, 24482,
                                                                       10140, 10200, 16820,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 25679, 0, 6,
                                                                       24482, 15605, 24545,
                                                                       10200, 10260, 16910,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsh_three_center_electron_repulsion_0(buffer, 25805, 0, 6,
                                                                       24545, 15650, 24608,
                                                                       10260, 10320, 17000,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 25931, 0, 3, 6,
                                                                       24671, 16010, 24860,
                                                                       25427, 16730, 25553,
                                                                       10440, 10620, 17630,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 26309, 0, 3, 6,
                                                                       24860, 16145, 25049,
                                                                       25553, 16820, 25679,
                                                                       10620, 10800, 17900,
                                                                       ncols, gamma, p, q);

                    compute_prim_dph_three_center_electron_repulsion_0(buffer, 26687, 0, 3, 6,
                                                                       25049, 16280, 25238,
                                                                       25679, 16910, 25805,
                                                                       10800, 10980, 18170,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 27065, 0, 6,
                                                                       25427, 16730, 25553,
                                                                       11340, 11440, 18740,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 27275, 0, 6,
                                                                       25553, 16820, 25679,
                                                                       11440, 11540, 18890,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsh_three_center_electron_repulsion_0(buffer, 27485, 0, 6,
                                                                       25679, 16910, 25805,
                                                                       11540, 11640, 19040,
                                                                       ncols, gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 27695, 0, 3, 6,
                                                                       24671, 24860, 25931,
                                                                       17630, 26309, 27065,
                                                                       18740, 27275, 11840,
                                                                       12140, 20090, ncols,
                                                                       gamma, p, q);

                    compute_prim_fph_three_center_electron_repulsion_0(buffer, 28325, 0, 3, 6,
                                                                       24860, 25049, 26309,
                                                                       17900, 26687, 27275,
                                                                       18890, 27485, 12140,
                                                                       12440, 20540, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 28955, 0, 6,
                                                                       27065, 18740, 27275,
                                                                       13040, 13190, 21440,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsh_three_center_electron_repulsion_0(buffer, 29270, 0, 6,
                                                                       27275, 18890, 27485,
                                                                       13190, 13340, 21665,
                                                                       ncols, gamma, p, q);

                    compute_prim_gph_three_center_electron_repulsion_0(buffer, 29585, 0, 3, 6,
                                                                       25931, 26309, 27695,
                                                                       20090, 28325, 28955,
                                                                       21440, 29270, 13640,
                                                                       14090, 23240, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 30530, 6, 14990,
                                                                       15005, 23915, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 30558, 6, 15005,
                                                                       15020, 23936, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 30586, 6, 15020,
                                                                       15035, 23957, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 30614, 6, 15035,
                                                                       15050, 23978, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 30642, 6, 15050,
                                                                       15065, 23999, ncols,
                                                                       gamma, p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 30670, 6, 15065,
                                                                       15080, 24020, ncols,
                                                                       gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 30698, 3, 6,
                                                                       30530, 23915, 30558,
                                                                       15110, 15155, 24041,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 30782, 3, 6,
                                                                       30558, 23936, 30586,
                                                                       15155, 15200, 24104,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 30866, 3, 6,
                                                                       30586, 23957, 30614,
                                                                       15200, 15245, 24167,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 30950, 3, 6,
                                                                       30614, 23978, 30642,
                                                                       15245, 15290, 24230,
                                                                       ncols, gamma, p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 31034, 3, 6,
                                                                       30642, 23999, 30670,
                                                                       15290, 15335, 24293,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 31118, 0, 6,
                                                                       30530, 23915, 30558,
                                                                       15425, 15470, 24356,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 31202, 0, 6,
                                                                       30558, 23936, 30586,
                                                                       15470, 15515, 24419,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 31286, 0, 6,
                                                                       30586, 23957, 30614,
                                                                       15515, 15560, 24482,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 31370, 0, 6,
                                                                       30614, 23978, 30642,
                                                                       15560, 15605, 24545,
                                                                       ncols, gamma, p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 31454, 0, 6,
                                                                       30642, 23999, 30670,
                                                                       15605, 15650, 24608,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 31538, 0, 3, 6,
                                                                       30698, 24041, 30782,
                                                                       31118, 24356, 31202,
                                                                       15740, 15875, 24671,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 31790, 0, 3, 6,
                                                                       30782, 24104, 30866,
                                                                       31202, 24419, 31286,
                                                                       15875, 16010, 24860,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 32042, 0, 3, 6,
                                                                       30866, 24167, 30950,
                                                                       31286, 24482, 31370,
                                                                       16010, 16145, 25049,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 32294, 0, 3, 6,
                                                                       30950, 24230, 31034,
                                                                       31370, 24545, 31454,
                                                                       16145, 16280, 25238,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 32546, 0, 6,
                                                                       31118, 24356, 31202,
                                                                       16550, 16640, 25427,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 32714, 0, 6,
                                                                       31202, 24419, 31286,
                                                                       16640, 16730, 25553,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 32882, 0, 6,
                                                                       31286, 24482, 31370,
                                                                       16730, 16820, 25679,
                                                                       ncols, gamma, p, q);

                    compute_prim_dsi_three_center_electron_repulsion_0(buffer, 33050, 0, 6,
                                                                       31370, 24545, 31454,
                                                                       16820, 16910, 25805,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 33218, 0, 3, 6,
                                                                       31538, 24671, 31790,
                                                                       32546, 25427, 32714,
                                                                       17090, 17360, 25931,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 33722, 0, 3, 6,
                                                                       31790, 24860, 32042,
                                                                       32714, 25553, 32882,
                                                                       17360, 17630, 26309,
                                                                       ncols, gamma, p, q);

                    compute_prim_dpi_three_center_electron_repulsion_0(buffer, 34226, 0, 3, 6,
                                                                       32042, 25049, 32294,
                                                                       32882, 25679, 33050,
                                                                       17630, 17900, 26687,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 34730, 0, 6,
                                                                       32546, 25427, 32714,
                                                                       18440, 18590, 27065,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 35010, 0, 6,
                                                                       32714, 25553, 32882,
                                                                       18590, 18740, 27275,
                                                                       ncols, gamma, p, q);

                    compute_prim_fsi_three_center_electron_repulsion_0(buffer, 35290, 0, 6,
                                                                       32882, 25679, 33050,
                                                                       18740, 18890, 27485,
                                                                       ncols, gamma, p, q);

                    compute_prim_fpi_three_center_electron_repulsion_0(buffer, 35570, 0, 3, 6,
                                                                       31538, 31790, 33218,
                                                                       25931, 33722, 34730,
                                                                       27065, 35010, 19190,
                                                                       19640, 27695, ncols,
                                                                       gamma, p, q);

                    compute_prim_fpi_three_center_electron_repulsion_0(buffer, 36410, 0, 3, 6,
                                                                       31790, 32042, 33722,
                                                                       26309, 34226, 35010,
                                                                       27275, 35290, 19640,
                                                                       20090, 28325, ncols,
                                                                       gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 37250, 0, 6,
                                                                       34730, 27065, 35010,
                                                                       20990, 21215, 28955,
                                                                       ncols, gamma, p, q);

                    compute_prim_gsi_three_center_electron_repulsion_0(buffer, 37670, 0, 6,
                                                                       35010, 27275, 35290,
                                                                       21215, 21440, 29270,
                                                                       ncols, gamma, p, q);

                    compute_prim_gpi_three_center_electron_repulsion_0(buffer, 38090, 0, 3, 6,
                                                                       33218, 33722, 35570,
                                                                       27695, 36410, 37250,
                                                                       28955, 37670, 21890,
                                                                       22565, 29585, ncols,
                                                                       gamma, p, q);

                    simdgeo::geom_s_x(buffer, 39350, 38090, 15, 28, ncols, beta);

                    simdgeo::geom_s_y(buffer, 39770, 38090, 15, 28, ncols, beta);

                    simdgeo::geom_s_z(buffer, 40190, 38090, 15, 28, ncols, beta);

                    simdfunc::contract_primitives(buffer, 40610, 39350, 1260, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 41870, 40610, 15, 1, nmax);

        simdtrf::transform_g_outer(values + n * npairs, nvalues, buffer, 41870, 13, nmax);

        simdtrf::transform_i_inner(buffer, 41870, 41030, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 117 * nvalues + n * npairs, nvalues, buffer, 41870,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 41870, 41450, 15, 1, nmax);

        simdtrf::transform_g_outer(values + 234 * nvalues + n * npairs, nvalues, buffer, 41870,
                                   13, nmax);
    }

    for (size_t m = 0; m < 351; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
