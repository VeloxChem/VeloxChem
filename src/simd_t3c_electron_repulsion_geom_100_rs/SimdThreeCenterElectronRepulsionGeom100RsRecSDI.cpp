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


#include "SimdThreeCenterElectronRepulsionGeom100RsRecSDI.hpp"

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
#include "SimdTransformD.hpp"
#include "SimdTransformI.hpp"

namespace simdt3ceri {  // simdt3ceri namespace

auto
compute_rs_geom_100_sdi_three_center_electron_repulsion(double               *values,
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
            false, std::string("compute_rs_geom_100_sdi_three_center_electron_repulsion: Number of values exceeds number of atom pairs"));
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

    const auto nmax = simdfunc::prepare_buffer(buffer, 21213, 0, 0, dimensions);

    if (nmax == 0)
    {
        std::fill(values, values + 390 * natoms * npairs, 0.0);

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
        simdfunc::prepare_buffer(buffer, 21213, 20127, 1008, dimensions);

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

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 181, 0, 6, 12, 13,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 184, 0, 6, 13, 14,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 187, 0, 6, 14, 15,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 190, 0, 6, 15, 16,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 193, 0, 6, 16, 17,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 196, 0, 6, 17, 18,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 199, 0, 6, 18, 19,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 202, 0, 6, 23, 24,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 205, 0, 6, 24, 25,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 208, 0, 6, 25, 26,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 211, 0, 6, 26, 27,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 214, 0, 6, 27, 28,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 217, 0, 6, 28, 29,
                                                                       ncols, gamma, q);

                    compute_prim_pss_three_center_electron_repulsion_0(buffer, 220, 0, 6, 29, 30,
                                                                       ncols, gamma, q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 223, 0, 6, 10, 11,
                                                                       31, 34, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 232, 0, 6, 11, 12,
                                                                       34, 37, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 241, 0, 6, 12, 13,
                                                                       37, 40, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 250, 0, 6, 13, 14,
                                                                       40, 43, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 259, 0, 6, 14, 15,
                                                                       43, 46, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 268, 0, 6, 15, 16,
                                                                       46, 49, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 277, 0, 6, 16, 17,
                                                                       49, 52, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 286, 0, 6, 17, 18,
                                                                       52, 55, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 295, 0, 6, 21, 22,
                                                                       58, 61, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 304, 0, 6, 22, 23,
                                                                       61, 64, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 313, 0, 6, 23, 24,
                                                                       64, 67, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 322, 0, 6, 24, 25,
                                                                       67, 70, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 331, 0, 6, 25, 26,
                                                                       70, 73, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 340, 0, 6, 26, 27,
                                                                       73, 76, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 349, 0, 6, 27, 28,
                                                                       76, 79, ncols, gamma, p,
                                                                       q);

                    compute_prim_pps_three_center_electron_repulsion_0(buffer, 358, 0, 6, 28, 29,
                                                                       79, 82, ncols, gamma, p,
                                                                       q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 367, 0, 3, 6, 31,
                                                                       34, 85, 91, 223, 232,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 385, 0, 3, 6, 34,
                                                                       37, 91, 97, 232, 241,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 403, 0, 3, 6, 37,
                                                                       40, 97, 103, 241, 250,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 421, 0, 3, 6, 40,
                                                                       43, 103, 109, 250, 259,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 439, 0, 3, 6, 43,
                                                                       46, 109, 115, 259, 268,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 457, 0, 3, 6, 46,
                                                                       49, 115, 121, 268, 277,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 475, 0, 3, 6, 49,
                                                                       52, 121, 127, 277, 286,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 493, 0, 3, 6, 58,
                                                                       61, 133, 139, 295, 304,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 511, 0, 3, 6, 61,
                                                                       64, 139, 145, 304, 313,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 529, 0, 3, 6, 64,
                                                                       67, 145, 151, 313, 322,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 547, 0, 3, 6, 67,
                                                                       70, 151, 157, 322, 331,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 565, 0, 3, 6, 70,
                                                                       73, 157, 163, 331, 340,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 583, 0, 3, 6, 73,
                                                                       76, 163, 169, 340, 349,
                                                                       ncols, gamma, p, q);

                    compute_prim_pds_three_center_electron_repulsion_0(buffer, 601, 0, 3, 6, 76,
                                                                       79, 169, 175, 349, 358,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 619, 6, 12, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 622, 6, 13, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 625, 6, 14, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 628, 6, 15, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 631, 6, 16, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 634, 6, 17, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 637, 6, 18, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 640, 6, 19, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 643, 6, 23, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 646, 6, 24, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 649, 6, 25, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 652, 6, 26, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 655, 6, 27, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 658, 6, 28, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 661, 6, 29, ncols,
                                                                       p, q);

                    compute_prim_ssp_three_center_electron_repulsion_0(buffer, 664, 6, 30, ncols,
                                                                       p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 667, 6, 12, 37,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 676, 6, 13, 40,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 685, 6, 14, 43,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 694, 6, 15, 46,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 703, 6, 16, 49,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 712, 6, 17, 52,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 721, 6, 18, 55,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 730, 6, 23, 64,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 739, 6, 24, 67,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 748, 6, 25, 70,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 757, 6, 26, 73,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 766, 6, 27, 76,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 775, 6, 28, 79,
                                                                       ncols, p, q);

                    compute_prim_spp_three_center_electron_repulsion_0(buffer, 784, 6, 29, 82,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 793, 6, 37, 97,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 811, 6, 40, 103,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 829, 6, 43, 109,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 847, 6, 46, 115,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 865, 6, 49, 121,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 883, 6, 52, 127,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 901, 6, 64, 145,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 919, 6, 67, 151,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 937, 6, 70, 157,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 955, 6, 73, 163,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 973, 6, 76, 169,
                                                                       ncols, p, q);

                    compute_prim_sdp_three_center_electron_repulsion_0(buffer, 991, 6, 79, 175,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1009, 6, 12, 181,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1018, 6, 13, 184,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1027, 6, 14, 187,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1036, 6, 15, 190,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1045, 6, 16, 193,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1054, 6, 17, 196,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1063, 6, 18, 199,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1072, 6, 23, 202,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1081, 6, 24, 205,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1090, 6, 25, 208,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1099, 6, 26, 211,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1108, 6, 27, 214,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1117, 6, 28, 217,
                                                                       ncols, p, q);

                    compute_prim_psp_three_center_electron_repulsion_0(buffer, 1126, 6, 29, 220,
                                                                       ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1135, 6, 37, 181,
                                                                       241, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1162, 6, 40, 184,
                                                                       250, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1189, 6, 43, 187,
                                                                       259, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1216, 6, 46, 190,
                                                                       268, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1243, 6, 49, 193,
                                                                       277, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1270, 6, 52, 196,
                                                                       286, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1297, 6, 64, 202,
                                                                       313, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1324, 6, 67, 205,
                                                                       322, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1351, 6, 70, 208,
                                                                       331, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1378, 6, 73, 211,
                                                                       340, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1405, 6, 76, 214,
                                                                       349, ncols, p, q);

                    compute_prim_ppp_three_center_electron_repulsion_0(buffer, 1432, 6, 79, 217,
                                                                       358, ncols, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1459, 3, 6, 97,
                                                                       1135, 241, 1162, 403,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1513, 3, 6, 103,
                                                                       1162, 250, 1189, 421,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1567, 3, 6, 109,
                                                                       1189, 259, 1216, 439,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1621, 3, 6, 115,
                                                                       1216, 268, 1243, 457,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1675, 3, 6, 121,
                                                                       1243, 277, 1270, 475,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1729, 3, 6, 145,
                                                                       1297, 313, 1324, 529,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1783, 3, 6, 151,
                                                                       1324, 322, 1351, 547,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1837, 3, 6, 157,
                                                                       1351, 331, 1378, 565,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1891, 3, 6, 163,
                                                                       1378, 340, 1405, 583,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdp_three_center_electron_repulsion_0(buffer, 1945, 3, 6, 169,
                                                                       1405, 349, 1432, 601,
                                                                       ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 1999, 6, 10, 11,
                                                                       619, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2005, 6, 11, 12,
                                                                       622, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2011, 6, 12, 13,
                                                                       625, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2017, 6, 13, 14,
                                                                       628, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2023, 6, 14, 15,
                                                                       631, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2029, 6, 15, 16,
                                                                       634, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2035, 6, 16, 17,
                                                                       637, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2041, 6, 17, 18,
                                                                       640, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2047, 6, 21, 22,
                                                                       643, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2053, 6, 22, 23,
                                                                       646, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2059, 6, 23, 24,
                                                                       649, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2065, 6, 24, 25,
                                                                       652, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2071, 6, 25, 26,
                                                                       655, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2077, 6, 26, 27,
                                                                       658, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2083, 6, 27, 28,
                                                                       661, ncols, gamma, p, q);

                    compute_prim_ssd_three_center_electron_repulsion_0(buffer, 2089, 6, 28, 29,
                                                                       664, ncols, gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2095, 3, 6, 1999,
                                                                       619, 2005, 667, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2113, 3, 6, 2005,
                                                                       622, 2011, 676, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2131, 3, 6, 2011,
                                                                       625, 2017, 685, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2149, 3, 6, 2017,
                                                                       628, 2023, 694, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2167, 3, 6, 2023,
                                                                       631, 2029, 703, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2185, 3, 6, 2029,
                                                                       634, 2035, 712, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2203, 3, 6, 2035,
                                                                       637, 2041, 721, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2221, 3, 6, 2047,
                                                                       643, 2053, 730, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2239, 3, 6, 2053,
                                                                       646, 2059, 739, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2257, 3, 6, 2059,
                                                                       649, 2065, 748, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2275, 3, 6, 2065,
                                                                       652, 2071, 757, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2293, 3, 6, 2071,
                                                                       655, 2077, 766, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2311, 3, 6, 2077,
                                                                       658, 2083, 775, ncols,
                                                                       gamma, p, q);

                    compute_prim_spd_three_center_electron_repulsion_0(buffer, 2329, 3, 6, 2083,
                                                                       661, 2089, 784, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2347, 3, 6, 2095,
                                                                       667, 2113, 85, 91, 793,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2383, 3, 6, 2113,
                                                                       676, 2131, 91, 97, 811,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2419, 3, 6, 2131,
                                                                       685, 2149, 97, 103, 829,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2455, 3, 6, 2149,
                                                                       694, 2167, 103, 109, 847,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2491, 3, 6, 2167,
                                                                       703, 2185, 109, 115, 865,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2527, 3, 6, 2185,
                                                                       712, 2203, 115, 121, 883,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2563, 3, 6, 2221,
                                                                       730, 2239, 133, 139, 901,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2599, 3, 6, 2239,
                                                                       739, 2257, 139, 145, 919,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2635, 3, 6, 2257,
                                                                       748, 2275, 145, 151, 937,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2671, 3, 6, 2275,
                                                                       757, 2293, 151, 157, 955,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2707, 3, 6, 2293,
                                                                       766, 2311, 157, 163, 973,
                                                                       ncols, gamma, p, q);

                    compute_prim_sdd_three_center_electron_repulsion_0(buffer, 2743, 3, 6, 2311,
                                                                       775, 2329, 163, 169, 991,
                                                                       ncols, gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2779, 0, 6, 1999,
                                                                       619, 2005, 1009, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2797, 0, 6, 2005,
                                                                       622, 2011, 1018, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2815, 0, 6, 2011,
                                                                       625, 2017, 1027, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2833, 0, 6, 2017,
                                                                       628, 2023, 1036, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2851, 0, 6, 2023,
                                                                       631, 2029, 1045, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2869, 0, 6, 2029,
                                                                       634, 2035, 1054, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2887, 0, 6, 2035,
                                                                       637, 2041, 1063, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2905, 0, 6, 2047,
                                                                       643, 2053, 1072, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2923, 0, 6, 2053,
                                                                       646, 2059, 1081, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2941, 0, 6, 2059,
                                                                       649, 2065, 1090, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2959, 0, 6, 2065,
                                                                       652, 2071, 1099, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2977, 0, 6, 2071,
                                                                       655, 2077, 1108, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 2995, 0, 6, 2077,
                                                                       658, 2083, 1117, ncols,
                                                                       gamma, p, q);

                    compute_prim_psd_three_center_electron_repulsion_0(buffer, 3013, 0, 6, 2083,
                                                                       661, 2089, 1126, ncols,
                                                                       gamma, p, q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3031, 0, 3, 6,
                                                                       2095, 667, 2113, 2779,
                                                                       1009, 2797, 223, 232,
                                                                       1135, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3085, 0, 3, 6,
                                                                       2113, 676, 2131, 2797,
                                                                       1018, 2815, 232, 241,
                                                                       1162, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3139, 0, 3, 6,
                                                                       2131, 685, 2149, 2815,
                                                                       1027, 2833, 241, 250,
                                                                       1189, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3193, 0, 3, 6,
                                                                       2149, 694, 2167, 2833,
                                                                       1036, 2851, 250, 259,
                                                                       1216, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3247, 0, 3, 6,
                                                                       2167, 703, 2185, 2851,
                                                                       1045, 2869, 259, 268,
                                                                       1243, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3301, 0, 3, 6,
                                                                       2185, 712, 2203, 2869,
                                                                       1054, 2887, 268, 277,
                                                                       1270, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3355, 0, 3, 6,
                                                                       2221, 730, 2239, 2905,
                                                                       1072, 2923, 295, 304,
                                                                       1297, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3409, 0, 3, 6,
                                                                       2239, 739, 2257, 2923,
                                                                       1081, 2941, 304, 313,
                                                                       1324, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3463, 0, 3, 6,
                                                                       2257, 748, 2275, 2941,
                                                                       1090, 2959, 313, 322,
                                                                       1351, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3517, 0, 3, 6,
                                                                       2275, 757, 2293, 2959,
                                                                       1099, 2977, 322, 331,
                                                                       1378, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3571, 0, 3, 6,
                                                                       2293, 766, 2311, 2977,
                                                                       1108, 2995, 331, 340,
                                                                       1405, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppd_three_center_electron_repulsion_0(buffer, 3625, 0, 3, 6,
                                                                       2311, 775, 2329, 2995,
                                                                       1117, 3013, 340, 349,
                                                                       1432, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3679, 0, 3, 6,
                                                                       2347, 793, 2383, 3031,
                                                                       1135, 3085, 367, 385,
                                                                       1459, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3787, 0, 3, 6,
                                                                       2383, 811, 2419, 3085,
                                                                       1162, 3139, 385, 403,
                                                                       1513, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 3895, 0, 3, 6,
                                                                       2419, 829, 2455, 3139,
                                                                       1189, 3193, 403, 421,
                                                                       1567, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4003, 0, 3, 6,
                                                                       2455, 847, 2491, 3193,
                                                                       1216, 3247, 421, 439,
                                                                       1621, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4111, 0, 3, 6,
                                                                       2491, 865, 2527, 3247,
                                                                       1243, 3301, 439, 457,
                                                                       1675, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4219, 0, 3, 6,
                                                                       2563, 901, 2599, 3355,
                                                                       1297, 3409, 493, 511,
                                                                       1729, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4327, 0, 3, 6,
                                                                       2599, 919, 2635, 3409,
                                                                       1324, 3463, 511, 529,
                                                                       1783, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4435, 0, 3, 6,
                                                                       2635, 937, 2671, 3463,
                                                                       1351, 3517, 529, 547,
                                                                       1837, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4543, 0, 3, 6,
                                                                       2671, 955, 2707, 3517,
                                                                       1378, 3571, 547, 565,
                                                                       1891, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdd_three_center_electron_repulsion_0(buffer, 4651, 0, 3, 6,
                                                                       2707, 973, 2743, 3571,
                                                                       1405, 3625, 565, 583,
                                                                       1945, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4759, 6, 619, 622,
                                                                       2011, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4769, 6, 622, 625,
                                                                       2017, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4779, 6, 625, 628,
                                                                       2023, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4789, 6, 628, 631,
                                                                       2029, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4799, 6, 631, 634,
                                                                       2035, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4809, 6, 634, 637,
                                                                       2041, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4819, 6, 643, 646,
                                                                       2059, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4829, 6, 646, 649,
                                                                       2065, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4839, 6, 649, 652,
                                                                       2071, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4849, 6, 652, 655,
                                                                       2077, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4859, 6, 655, 658,
                                                                       2083, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssf_three_center_electron_repulsion_0(buffer, 4869, 6, 658, 661,
                                                                       2089, ncols, gamma, p,
                                                                       q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4879, 3, 6, 4759,
                                                                       2011, 4769, 2131, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4909, 3, 6, 4769,
                                                                       2017, 4779, 2149, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4939, 3, 6, 4779,
                                                                       2023, 4789, 2167, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4969, 3, 6, 4789,
                                                                       2029, 4799, 2185, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 4999, 3, 6, 4799,
                                                                       2035, 4809, 2203, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5029, 3, 6, 4819,
                                                                       2059, 4829, 2257, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5059, 3, 6, 4829,
                                                                       2065, 4839, 2275, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5089, 3, 6, 4839,
                                                                       2071, 4849, 2293, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5119, 3, 6, 4849,
                                                                       2077, 4859, 2311, ncols,
                                                                       gamma, p, q);

                    compute_prim_spf_three_center_electron_repulsion_0(buffer, 5149, 3, 6, 4859,
                                                                       2083, 4869, 2329, ncols,
                                                                       gamma, p, q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5179, 3, 6, 4879,
                                                                       2131, 4909, 793, 811,
                                                                       2419, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5239, 3, 6, 4909,
                                                                       2149, 4939, 811, 829,
                                                                       2455, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5299, 3, 6, 4939,
                                                                       2167, 4969, 829, 847,
                                                                       2491, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5359, 3, 6, 4969,
                                                                       2185, 4999, 847, 865,
                                                                       2527, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5419, 3, 6, 5029,
                                                                       2257, 5059, 901, 919,
                                                                       2635, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5479, 3, 6, 5059,
                                                                       2275, 5089, 919, 937,
                                                                       2671, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5539, 3, 6, 5089,
                                                                       2293, 5119, 937, 955,
                                                                       2707, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdf_three_center_electron_repulsion_0(buffer, 5599, 3, 6, 5119,
                                                                       2311, 5149, 955, 973,
                                                                       2743, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5659, 0, 6, 4759,
                                                                       2011, 4769, 1009, 1018,
                                                                       2815, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5689, 0, 6, 4769,
                                                                       2017, 4779, 1018, 1027,
                                                                       2833, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5719, 0, 6, 4779,
                                                                       2023, 4789, 1027, 1036,
                                                                       2851, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5749, 0, 6, 4789,
                                                                       2029, 4799, 1036, 1045,
                                                                       2869, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5779, 0, 6, 4799,
                                                                       2035, 4809, 1045, 1054,
                                                                       2887, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5809, 0, 6, 4819,
                                                                       2059, 4829, 1072, 1081,
                                                                       2941, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5839, 0, 6, 4829,
                                                                       2065, 4839, 1081, 1090,
                                                                       2959, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5869, 0, 6, 4839,
                                                                       2071, 4849, 1090, 1099,
                                                                       2977, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5899, 0, 6, 4849,
                                                                       2077, 4859, 1099, 1108,
                                                                       2995, ncols, gamma, p,
                                                                       q);

                    compute_prim_psf_three_center_electron_repulsion_0(buffer, 5929, 0, 6, 4859,
                                                                       2083, 4869, 1108, 1117,
                                                                       3013, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 5959, 0, 3, 6,
                                                                       4879, 2131, 4909, 5659,
                                                                       2815, 5689, 1135, 1162,
                                                                       3139, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6049, 0, 3, 6,
                                                                       4909, 2149, 4939, 5689,
                                                                       2833, 5719, 1162, 1189,
                                                                       3193, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6139, 0, 3, 6,
                                                                       4939, 2167, 4969, 5719,
                                                                       2851, 5749, 1189, 1216,
                                                                       3247, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6229, 0, 3, 6,
                                                                       4969, 2185, 4999, 5749,
                                                                       2869, 5779, 1216, 1243,
                                                                       3301, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6319, 0, 3, 6,
                                                                       5029, 2257, 5059, 5809,
                                                                       2941, 5839, 1297, 1324,
                                                                       3463, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6409, 0, 3, 6,
                                                                       5059, 2275, 5089, 5839,
                                                                       2959, 5869, 1324, 1351,
                                                                       3517, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6499, 0, 3, 6,
                                                                       5089, 2293, 5119, 5869,
                                                                       2977, 5899, 1351, 1378,
                                                                       3571, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppf_three_center_electron_repulsion_0(buffer, 6589, 0, 3, 6,
                                                                       5119, 2311, 5149, 5899,
                                                                       2995, 5929, 1378, 1405,
                                                                       3625, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 6679, 0, 3, 6,
                                                                       5179, 2419, 5239, 5959,
                                                                       3139, 6049, 1459, 1513,
                                                                       3895, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 6859, 0, 3, 6,
                                                                       5239, 2455, 5299, 6049,
                                                                       3193, 6139, 1513, 1567,
                                                                       4003, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 7039, 0, 3, 6,
                                                                       5299, 2491, 5359, 6139,
                                                                       3247, 6229, 1567, 1621,
                                                                       4111, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 7219, 0, 3, 6,
                                                                       5419, 2635, 5479, 6319,
                                                                       3463, 6409, 1729, 1783,
                                                                       4435, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 7399, 0, 3, 6,
                                                                       5479, 2671, 5539, 6409,
                                                                       3517, 6499, 1783, 1837,
                                                                       4543, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdf_three_center_electron_repulsion_0(buffer, 7579, 0, 3, 6,
                                                                       5539, 2707, 5599, 6499,
                                                                       3571, 6589, 1837, 1891,
                                                                       4651, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7759, 6, 1999,
                                                                       2005, 4759, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7774, 6, 2005,
                                                                       2011, 4769, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7789, 6, 2011,
                                                                       2017, 4779, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7804, 6, 2017,
                                                                       2023, 4789, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7819, 6, 2023,
                                                                       2029, 4799, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7834, 6, 2029,
                                                                       2035, 4809, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7849, 6, 2047,
                                                                       2053, 4819, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7864, 6, 2053,
                                                                       2059, 4829, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7879, 6, 2059,
                                                                       2065, 4839, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7894, 6, 2065,
                                                                       2071, 4849, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7909, 6, 2071,
                                                                       2077, 4859, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssg_three_center_electron_repulsion_0(buffer, 7924, 6, 2077,
                                                                       2083, 4869, ncols, gamma,
                                                                       p, q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7939, 3, 6, 7759,
                                                                       4759, 7774, 2095, 2113,
                                                                       4879, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 7984, 3, 6, 7774,
                                                                       4769, 7789, 2113, 2131,
                                                                       4909, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8029, 3, 6, 7789,
                                                                       4779, 7804, 2131, 2149,
                                                                       4939, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8074, 3, 6, 7804,
                                                                       4789, 7819, 2149, 2167,
                                                                       4969, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8119, 3, 6, 7819,
                                                                       4799, 7834, 2167, 2185,
                                                                       4999, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8164, 3, 6, 7849,
                                                                       4819, 7864, 2221, 2239,
                                                                       5029, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8209, 3, 6, 7864,
                                                                       4829, 7879, 2239, 2257,
                                                                       5059, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8254, 3, 6, 7879,
                                                                       4839, 7894, 2257, 2275,
                                                                       5089, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8299, 3, 6, 7894,
                                                                       4849, 7909, 2275, 2293,
                                                                       5119, ncols, gamma, p,
                                                                       q);

                    compute_prim_spg_three_center_electron_repulsion_0(buffer, 8344, 3, 6, 7909,
                                                                       4859, 7924, 2293, 2311,
                                                                       5149, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8389, 3, 6, 7939,
                                                                       4879, 7984, 2347, 2383,
                                                                       5179, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8479, 3, 6, 7984,
                                                                       4909, 8029, 2383, 2419,
                                                                       5239, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8569, 3, 6, 8029,
                                                                       4939, 8074, 2419, 2455,
                                                                       5299, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8659, 3, 6, 8074,
                                                                       4969, 8119, 2455, 2491,
                                                                       5359, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8749, 3, 6, 8164,
                                                                       5029, 8209, 2563, 2599,
                                                                       5419, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8839, 3, 6, 8209,
                                                                       5059, 8254, 2599, 2635,
                                                                       5479, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 8929, 3, 6, 8254,
                                                                       5089, 8299, 2635, 2671,
                                                                       5539, ncols, gamma, p,
                                                                       q);

                    compute_prim_sdg_three_center_electron_repulsion_0(buffer, 9019, 3, 6, 8299,
                                                                       5119, 8344, 2671, 2707,
                                                                       5599, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9109, 0, 6, 7759,
                                                                       4759, 7774, 2779, 2797,
                                                                       5659, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9154, 0, 6, 7774,
                                                                       4769, 7789, 2797, 2815,
                                                                       5689, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9199, 0, 6, 7789,
                                                                       4779, 7804, 2815, 2833,
                                                                       5719, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9244, 0, 6, 7804,
                                                                       4789, 7819, 2833, 2851,
                                                                       5749, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9289, 0, 6, 7819,
                                                                       4799, 7834, 2851, 2869,
                                                                       5779, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9334, 0, 6, 7849,
                                                                       4819, 7864, 2905, 2923,
                                                                       5809, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9379, 0, 6, 7864,
                                                                       4829, 7879, 2923, 2941,
                                                                       5839, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9424, 0, 6, 7879,
                                                                       4839, 7894, 2941, 2959,
                                                                       5869, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9469, 0, 6, 7894,
                                                                       4849, 7909, 2959, 2977,
                                                                       5899, ncols, gamma, p,
                                                                       q);

                    compute_prim_psg_three_center_electron_repulsion_0(buffer, 9514, 0, 6, 7909,
                                                                       4859, 7924, 2977, 2995,
                                                                       5929, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9559, 0, 3, 6,
                                                                       7939, 4879, 7984, 9109,
                                                                       5659, 9154, 3031, 3085,
                                                                       5959, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9694, 0, 3, 6,
                                                                       7984, 4909, 8029, 9154,
                                                                       5689, 9199, 3085, 3139,
                                                                       6049, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9829, 0, 3, 6,
                                                                       8029, 4939, 8074, 9199,
                                                                       5719, 9244, 3139, 3193,
                                                                       6139, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 9964, 0, 3, 6,
                                                                       8074, 4969, 8119, 9244,
                                                                       5749, 9289, 3193, 3247,
                                                                       6229, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 10099, 0, 3, 6,
                                                                       8164, 5029, 8209, 9334,
                                                                       5809, 9379, 3355, 3409,
                                                                       6319, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 10234, 0, 3, 6,
                                                                       8209, 5059, 8254, 9379,
                                                                       5839, 9424, 3409, 3463,
                                                                       6409, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 10369, 0, 3, 6,
                                                                       8254, 5089, 8299, 9424,
                                                                       5869, 9469, 3463, 3517,
                                                                       6499, ncols, gamma, p,
                                                                       q);

                    compute_prim_ppg_three_center_electron_repulsion_0(buffer, 10504, 0, 3, 6,
                                                                       8299, 5119, 8344, 9469,
                                                                       5899, 9514, 3517, 3571,
                                                                       6589, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 10639, 0, 3, 6,
                                                                       8389, 5179, 8479, 9559,
                                                                       5959, 9694, 3679, 3787,
                                                                       6679, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 10909, 0, 3, 6,
                                                                       8479, 5239, 8569, 9694,
                                                                       6049, 9829, 3787, 3895,
                                                                       6859, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 11179, 0, 3, 6,
                                                                       8569, 5299, 8659, 9829,
                                                                       6139, 9964, 3895, 4003,
                                                                       7039, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 11449, 0, 3, 6,
                                                                       8749, 5419, 8839, 10099,
                                                                       6319, 10234, 4219, 4327,
                                                                       7219, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 11719, 0, 3, 6,
                                                                       8839, 5479, 8929, 10234,
                                                                       6409, 10369, 4327, 4435,
                                                                       7399, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdg_three_center_electron_repulsion_0(buffer, 11989, 0, 3, 6,
                                                                       8929, 5539, 9019, 10369,
                                                                       6499, 10504, 4435, 4543,
                                                                       7579, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12259, 6, 4759,
                                                                       4769, 7789, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12280, 6, 4769,
                                                                       4779, 7804, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12301, 6, 4779,
                                                                       4789, 7819, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12322, 6, 4789,
                                                                       4799, 7834, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12343, 6, 4819,
                                                                       4829, 7879, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12364, 6, 4829,
                                                                       4839, 7894, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12385, 6, 4839,
                                                                       4849, 7909, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssh_three_center_electron_repulsion_0(buffer, 12406, 6, 4849,
                                                                       4859, 7924, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 12427, 3, 6,
                                                                       12259, 7789, 12280, 4879,
                                                                       4909, 8029, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 12490, 3, 6,
                                                                       12280, 7804, 12301, 4909,
                                                                       4939, 8074, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 12553, 3, 6,
                                                                       12301, 7819, 12322, 4939,
                                                                       4969, 8119, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 12616, 3, 6,
                                                                       12343, 7879, 12364, 5029,
                                                                       5059, 8254, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 12679, 3, 6,
                                                                       12364, 7894, 12385, 5059,
                                                                       5089, 8299, ncols, gamma,
                                                                       p, q);

                    compute_prim_sph_three_center_electron_repulsion_0(buffer, 12742, 3, 6,
                                                                       12385, 7909, 12406, 5089,
                                                                       5119, 8344, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 12805, 3, 6,
                                                                       12427, 8029, 12490, 5179,
                                                                       5239, 8569, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 12931, 3, 6,
                                                                       12490, 8074, 12553, 5239,
                                                                       5299, 8659, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 13057, 3, 6,
                                                                       12616, 8254, 12679, 5419,
                                                                       5479, 8929, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdh_three_center_electron_repulsion_0(buffer, 13183, 3, 6,
                                                                       12679, 8299, 12742, 5479,
                                                                       5539, 9019, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13309, 0, 6,
                                                                       12259, 7789, 12280, 5659,
                                                                       5689, 9199, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13372, 0, 6,
                                                                       12280, 7804, 12301, 5689,
                                                                       5719, 9244, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13435, 0, 6,
                                                                       12301, 7819, 12322, 5719,
                                                                       5749, 9289, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13498, 0, 6,
                                                                       12343, 7879, 12364, 5809,
                                                                       5839, 9424, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13561, 0, 6,
                                                                       12364, 7894, 12385, 5839,
                                                                       5869, 9469, ncols, gamma,
                                                                       p, q);

                    compute_prim_psh_three_center_electron_repulsion_0(buffer, 13624, 0, 6,
                                                                       12385, 7909, 12406, 5869,
                                                                       5899, 9514, ncols, gamma,
                                                                       p, q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 13687, 0, 3, 6,
                                                                       12427, 8029, 12490, 13309,
                                                                       9199, 13372, 5959, 6049,
                                                                       9829, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 13876, 0, 3, 6,
                                                                       12490, 8074, 12553, 13372,
                                                                       9244, 13435, 6049, 6139,
                                                                       9964, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 14065, 0, 3, 6,
                                                                       12616, 8254, 12679, 13498,
                                                                       9424, 13561, 6319, 6409,
                                                                       10369, ncols, gamma, p,
                                                                       q);

                    compute_prim_pph_three_center_electron_repulsion_0(buffer, 14254, 0, 3, 6,
                                                                       12679, 8299, 12742, 13561,
                                                                       9469, 13624, 6409, 6499,
                                                                       10504, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 14443, 0, 3, 6,
                                                                       12805, 8569, 12931, 13687,
                                                                       9829, 13876, 6679, 6859,
                                                                       11179, ncols, gamma, p,
                                                                       q);

                    compute_prim_pdh_three_center_electron_repulsion_0(buffer, 14821, 0, 3, 6,
                                                                       13057, 8929, 13183, 14065,
                                                                       10369, 14254, 7219, 7399,
                                                                       11989, ncols, gamma, p,
                                                                       q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15199, 6, 7759,
                                                                       7774, 12259, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15227, 6, 7774,
                                                                       7789, 12280, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15255, 6, 7789,
                                                                       7804, 12301, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15283, 6, 7804,
                                                                       7819, 12322, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15311, 6, 7849,
                                                                       7864, 12343, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15339, 6, 7864,
                                                                       7879, 12364, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15367, 6, 7879,
                                                                       7894, 12385, ncols, gamma,
                                                                       p, q);

                    compute_prim_ssi_three_center_electron_repulsion_0(buffer, 15395, 6, 7894,
                                                                       7909, 12406, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 15423, 3, 6,
                                                                       15199, 12259, 15227, 7939,
                                                                       7984, 12427, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 15507, 3, 6,
                                                                       15227, 12280, 15255, 7984,
                                                                       8029, 12490, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 15591, 3, 6,
                                                                       15255, 12301, 15283, 8029,
                                                                       8074, 12553, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 15675, 3, 6,
                                                                       15311, 12343, 15339, 8164,
                                                                       8209, 12616, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 15759, 3, 6,
                                                                       15339, 12364, 15367, 8209,
                                                                       8254, 12679, ncols, gamma,
                                                                       p, q);

                    compute_prim_spi_three_center_electron_repulsion_0(buffer, 15843, 3, 6,
                                                                       15367, 12385, 15395, 8254,
                                                                       8299, 12742, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 15927, 3, 6,
                                                                       15423, 12427, 15507, 8389,
                                                                       8479, 12805, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 16095, 3, 6,
                                                                       15507, 12490, 15591, 8479,
                                                                       8569, 12931, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 16263, 3, 6,
                                                                       15675, 12616, 15759, 8749,
                                                                       8839, 13057, ncols, gamma,
                                                                       p, q);

                    compute_prim_sdi_three_center_electron_repulsion_0(buffer, 16431, 3, 6,
                                                                       15759, 12679, 15843, 8839,
                                                                       8929, 13183, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16599, 0, 6,
                                                                       15199, 12259, 15227, 9109,
                                                                       9154, 13309, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16683, 0, 6,
                                                                       15227, 12280, 15255, 9154,
                                                                       9199, 13372, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16767, 0, 6,
                                                                       15255, 12301, 15283, 9199,
                                                                       9244, 13435, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16851, 0, 6,
                                                                       15311, 12343, 15339, 9334,
                                                                       9379, 13498, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 16935, 0, 6,
                                                                       15339, 12364, 15367, 9379,
                                                                       9424, 13561, ncols, gamma,
                                                                       p, q);

                    compute_prim_psi_three_center_electron_repulsion_0(buffer, 17019, 0, 6,
                                                                       15367, 12385, 15395, 9424,
                                                                       9469, 13624, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 17103, 0, 3, 6,
                                                                       15423, 12427, 15507,
                                                                       16599, 13309, 16683, 9559,
                                                                       9694, 13687, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 17355, 0, 3, 6,
                                                                       15507, 12490, 15591,
                                                                       16683, 13372, 16767, 9694,
                                                                       9829, 13876, ncols, gamma,
                                                                       p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 17607, 0, 3, 6,
                                                                       15675, 12616, 15759,
                                                                       16851, 13498, 16935,
                                                                       10099, 10234, 14065,
                                                                       ncols, gamma, p, q);

                    compute_prim_ppi_three_center_electron_repulsion_0(buffer, 17859, 0, 3, 6,
                                                                       15759, 12679, 15843,
                                                                       16935, 13561, 17019,
                                                                       10234, 10369, 14254,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 18111, 0, 3, 6,
                                                                       15927, 12805, 16095,
                                                                       17103, 13687, 17355,
                                                                       10639, 10909, 14443,
                                                                       ncols, gamma, p, q);

                    compute_prim_pdi_three_center_electron_repulsion_0(buffer, 18615, 0, 3, 6,
                                                                       16263, 13057, 16431,
                                                                       17607, 14065, 17859,
                                                                       11449, 11719, 14821,
                                                                       ncols, gamma, p, q);

                    simdgeo::geom_s_x(buffer, 19119, 18615, 1, 168, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 19287, 18615, 1, 168, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 19455, 18615, 1, 168, ncols, alpha);

                    simdgeo::geom_s_x(buffer, 19623, 18111, 1, 168, ncols, alpha);

                    simdgeo::geom_s_y(buffer, 19791, 18111, 1, 168, ncols, alpha);

                    simdgeo::geom_s_z(buffer, 19959, 18111, 1, 168, ncols, alpha);

                    simdfunc::contract_primitives(buffer, 20127, 19119, 1008, ncols);
                }
            }
        }

        simdtrf::transform_i_inner(buffer, 21135, 20127, 6, 1, nmax);

        simdtrf::transform_d_outer(values + n * npairs, nvalues, buffer, 21135, 13, nmax);

        simdtrf::transform_i_inner(buffer, 21135, 20295, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 65 * nvalues + n * npairs, nvalues, buffer, 21135,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 21135, 20463, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 130 * nvalues + n * npairs, nvalues, buffer, 21135,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 21135, 20631, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 195 * nvalues + n * npairs, nvalues, buffer, 21135,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 21135, 20799, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 260 * nvalues + n * npairs, nvalues, buffer, 21135,
                                   13, nmax);

        simdtrf::transform_i_inner(buffer, 21135, 20967, 6, 1, nmax);

        simdtrf::transform_d_outer(values + 325 * nvalues + n * npairs, nvalues, buffer, 21135,
                                   13, nmax);
    }

    for (size_t m = 0; m < 390; m++)
    {
        for (size_t n = 0; n < natoms; n++)
        {
            auto *pv = values + m * nvalues + n * npairs;

            std::fill(pv + nmax, pv + npairs, 0.0);
        }
    }
}

}  // namespace simdt3ceri
